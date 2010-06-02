/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Collapse short edges and combines edges that are in line.

    - collapse short edges. Length of edges to collapse provided as argument.
    - merge two edges if they are in line. Maximum angle provided as argument.
    - remove unused points.

    Cannot remove cells. Can remove faces and points but does not check
    for nonsense resulting topology.

    When collapsing an edge with one point on the boundary it will leave
    the boundary point intact. When both points inside it chooses random. When
    both points on boundary random again. Note: it should in fact use features
    where if one point is on a feature it collapses to that one. Alas we don't
    have features on a polyMesh.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "edgeCollapser.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "mathematicalConstants.H"
#include "PackedBoolList.H"
#include "SortableList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Get faceEdges in order of face points, i.e. faceEdges[0] is between
// f[0] and f[1]
labelList getSortedEdges
(
    const edgeList& edges,
    const labelList& f,
    const labelList& edgeLabels
)
{
    labelList faceEdges(edgeLabels.size(), -1);

    // Find starting pos in f for every edgeLabels
    forAll(edgeLabels, i)
    {
        label edgeI = edgeLabels[i];

        const edge& e = edges[edgeI];

        label fp = findIndex(f, e[0]);
        label fp1 = f.fcIndex(fp);

        if (f[fp1] == e[1])
        {
            // EdgeI between fp -> fp1
            faceEdges[fp] = edgeI;
        }
        else
        {
            // EdgeI between fp-1 -> fp
            faceEdges[f.rcIndex(fp)] = edgeI;
        }
    }

    return faceEdges;
}


// Merges edges which are in straight line. I.e. edge split by point.
label mergeEdges
(
    const polyMesh& mesh,
    const scalar maxCos,
    edgeCollapser& collapser
)
{
    const pointField& points = mesh.points();
    const edgeList& edges = mesh.edges();
    const labelListList& pointEdges = mesh.pointEdges();
    const labelList& region = collapser.pointRegion();
    const labelList& master = collapser.pointRegionMaster();

    label nCollapsed = 0;

    forAll(pointEdges, pointI)
    {
        const labelList& pEdges = pointEdges[pointI];

        if (pEdges.size() == 2)
        {
            const edge& leftE = edges[pEdges[0]];
            const edge& rightE = edges[pEdges[1]];

            // Get the two vertices on both sides of the point
            label leftV = leftE.otherVertex(pointI);
            label rightV = rightE.otherVertex(pointI);

            // Collapse only if none of the points part of merge network
            // or all of networks with different masters.
            label midMaster = -1;
            if (region[pointI] != -1)
            {
                midMaster = master[region[pointI]];
            }

            label leftMaster = -2;
            if (region[leftV] != -1)
            {
                leftMaster = master[region[leftV]];
            }

            label rightMaster = -3;
            if (region[rightV] != -1)
            {
                rightMaster = master[region[rightV]];
            }

            if
            (
                midMaster != leftMaster
             && midMaster != rightMaster
             && leftMaster != rightMaster
            )
            {
                // Check if the two edge are in line
                vector leftVec = points[pointI] - points[leftV];
                leftVec /= mag(leftVec) + VSMALL;

                vector rightVec = points[rightV] - points[pointI];
                rightVec /= mag(rightVec) + VSMALL;

                if ((leftVec & rightVec) > maxCos)
                {
                    // Collapse one (left) side of the edge. Make left vertex
                    // the master.
                    //if (collapser.unaffectedEdge(pEdges[0]))
                    {
                        collapser.collapseEdge(pEdges[0], leftV);
                        nCollapsed++;
                    }
                }
            }
        }
    }

    return nCollapsed;
}


// Return master point edge needs to be collapsed to (or -1)
label edgeMaster(const PackedBoolList& boundaryPoint, const edge& e)
{
    label masterPoint = -1;

    // Collapse edge to boundary point.
    if (boundaryPoint.get(e[0]))
    {
        if (boundaryPoint.get(e[1]))
        {
            // Both points on boundary. Choose one to collapse to.
            // Note: should look at feature edges/points!
            masterPoint = e[0];
        }
        else
        {
            masterPoint = e[0];
        }
    }
    else
    {
        if (boundaryPoint.get(e[1]))
        {
            masterPoint = e[1];
        }
        else
        {
            // None on boundary. Choose arbitrary.
            // Note: should look at geometry?
            masterPoint = e[0];
        }
    }
    return masterPoint;
}


label collapseSmallEdges
(
    const polyMesh& mesh,
    const PackedBoolList& boundaryPoint,
    const scalar minLen,
    edgeCollapser& collapser
)
{
    const pointField& points = mesh.points();
    const edgeList& edges = mesh.edges();

    // Collapse all edges that are too small. Choose intelligently which
    // point to collapse edge to.

    label nCollapsed = 0;

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        if (e.mag(points) < minLen)
        {
            label master = edgeMaster(boundaryPoint, e);

            if (master != -1) // && collapser.unaffectedEdge(edgeI))
            {
                collapser.collapseEdge(edgeI, master);
                nCollapsed++;
            }
        }
    }
    return nCollapsed;
}


// Faces which have edges just larger than collapse length but faces which
// are very small. This one tries to collapse them if it can be done with
// edge collapse. For faces where a face gets replace by two edges use
// collapseFaces
label collapseHighAspectFaces
(
    const polyMesh& mesh,
    const PackedBoolList& boundaryPoint,
    const scalar areaFac,
    const scalar edgeRatio,
    edgeCollapser& collapser
)
{
    const pointField& points = mesh.points();
    const edgeList& edges = mesh.edges();
    const faceList& faces = mesh.faces();
    const labelListList& faceEdges = mesh.faceEdges();

    scalarField magArea(mag(mesh.faceAreas()));

    label maxIndex = findMax(magArea);

    scalar minArea = areaFac * magArea[maxIndex];

    Info<< "Max face area:" << magArea[maxIndex] << endl
        << "Collapse area factor:" << areaFac << endl
        << "Collapse area:" << minArea << endl;

    label nCollapsed = 0;

    forAll(faces, faceI)
    {
        if (magArea[faceI] < minArea)
        {
            const face& f = faces[faceI];

            // Get the edges in face point order
            labelList fEdges(getSortedEdges(edges, f, faceEdges[faceI]));

            SortableList<scalar> lengths(fEdges.size());
            forAll(fEdges, i)
            {
                lengths[i] = edges[fEdges[i]].mag(points);
            }
            lengths.sort();


            label edgeI = -1;

            if (f.size() == 4)
            {
                // Compare second largest to smallest
                if (lengths[2] > edgeRatio*lengths[0])
                {
                    // Collapse smallest only. Triangle should be cleared
                    // next time around.
                    edgeI = fEdges[lengths.indices()[0]];
                }
            }
            else if (f.size() == 3)
            {
                // Compare second largest to smallest
                if (lengths[1] > edgeRatio*lengths[0])
                {
                    edgeI = fEdges[lengths.indices()[0]];
                }
            }


            if (edgeI != -1)
            {
                label master = edgeMaster(boundaryPoint, edges[edgeI]);

                if (master != -1)// && collapser.unaffectedEdge(edgeI))
                {
                    collapser.collapseEdge(edgeI, master);
                    nCollapsed++;
                }
            }
        }
    }

    return nCollapsed;
}


void set(const labelList& elems, const bool val, boolList& status)
{
    forAll(elems, i)
    {
        status[elems[i]] = val;
    }
}


// Tries to simplify polygons to face of minSize (4=quad, 3=triangle)
label simplifyFaces
(
    const polyMesh& mesh,
    const PackedBoolList& boundaryPoint,
    const label minSize,
    const scalar lenGap,
    edgeCollapser& collapser
)
{
    const pointField& points = mesh.points();
    const edgeList& edges = mesh.edges();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    const labelListList& faceEdges = mesh.faceEdges();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();
    const labelListList& pointCells = mesh.pointCells();
    const labelListList& cellEdges = mesh.cellEdges();

    label nCollapsed = 0;

    boolList protectedEdge(mesh.nEdges(), false);

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        if
        (
            f.size() > minSize
         && cells[faceOwner[faceI]].size() >= 6
         && (
                mesh.isInternalFace(faceI)
             && cells[faceNeighbour[faceI]].size() >= 6
            )
        )
        {
            // Get the edges in face point order
            labelList fEdges(getSortedEdges(edges, f, faceEdges[faceI]));

            SortableList<scalar> lengths(fEdges.size());
            forAll(fEdges, i)
            {
                lengths[i] = edges[fEdges[i]].mag(points);
            }
            lengths.sort();


            // Now find a gap in length between consecutive elements greater
            // than lenGap.

            label gapPos = -1;

            for (label i = f.size()-1-minSize; i >= 0; --i)
            {
                if (lengths[i+1] > lenGap*lengths[i])
                {
                    gapPos = i;

                    break;
                }
            }

            if (gapPos != -1)
            {
                //for (label i = gapPos; i >= 0; --i)
                label i = 0;  // Hack: collapse smallest edge only.
                {
                    label edgeI = fEdges[lengths.indices()[i]];

                    if (!protectedEdge[edgeI])
                    {
                        const edge& e = edges[edgeI];

                        label master = edgeMaster(boundaryPoint, e);

                        if (master != -1)
                        {
                            collapser.collapseEdge(edgeI, master);

                            // Protect all other edges on all cells using edge
                            // points.

                            const labelList& pCells0 = pointCells[e[0]];

                            forAll(pCells0, i)
                            {
                                set(cellEdges[pCells0[i]], true, protectedEdge);
                            }
                            const labelList& pCells1 = pointCells[e[1]];

                            forAll(pCells1, i)
                            {
                                set(cellEdges[pCells1[i]], true, protectedEdge);
                            }

                            nCollapsed++;
                        }
                    }
                }
            }
        }
    }

    return nCollapsed;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("overwrite", "");
    argList::validArgs.append("edge length [m]");
    argList::validArgs.append("merge angle (degrees)");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    scalar minLen(readScalar(IStringStream(args.additionalArgs()[0])()));
    scalar angle(readScalar(IStringStream(args.additionalArgs()[1])()));
    bool overwrite = args.optionFound("overwrite");

    scalar maxCos = Foam::cos(angle*mathematicalConstant::pi/180.0);

    Info<< "Merging:" << nl
        << "    edges with length less than " << minLen << " meters" << nl
        << "    edges split by a point with edges in line to within " << angle
        << " degrees" << nl
        << endl;


    bool meshChanged = false;

    while (true)
    {
        const faceList& faces = mesh.faces();

        // Get all points on the boundary
        PackedBoolList boundaryPoint(mesh.nPoints());

        label nIntFaces = mesh.nInternalFaces();
        for (label faceI = nIntFaces; faceI < mesh.nFaces(); faceI++)
        {
            const face& f = faces[faceI];

            forAll(f, fp)
            {
                boundaryPoint.set(f[fp], 1);
            }
        }

        // Edge collapsing engine
        edgeCollapser collapser(mesh);


        // Collapse all edges that are too small.
        label nCollapsed =
            collapseSmallEdges
            (
                mesh,
                boundaryPoint,
                minLen,
                collapser
            );
        Info<< "Collapsing " << nCollapsed << " small edges" << endl;


        // Remove midpoints on straight edges.
        if (nCollapsed == 0)
        {
            nCollapsed = mergeEdges(mesh, maxCos, collapser);
            Info<< "Collapsing " << nCollapsed << " in line edges" << endl;
        }


        // Remove small sliver faces that can be collapsed to single edge
        if (nCollapsed == 0)
        {
            nCollapsed =
                collapseHighAspectFaces
                (
                    mesh,
                    boundaryPoint,
                    1E-9,       // factor of largest face area
                    5,          // factor between smallest and largest edge on
                                // face
                    collapser
                );
            Info<< "Collapsing " << nCollapsed
                << " small high aspect ratio faces" << endl;
        }

        // Simplify faces to quads wherever possible
        //if (nCollapsed == 0)
        //{
        //    nCollapsed =
        //        simplifyFaces
        //        (
        //            mesh,
        //            boundaryPoint,
        //            4,              // minimum size of face
        //            0.2,            // gap in edge lengths on face
        //            collapser
        //        );
        //    Info<< "Collapsing " << nCollapsed << " polygonal faces" << endl;
        //}


        if (nCollapsed == 0)
        {
            break;
        }

        polyTopoChange meshMod(mesh);

        // Insert mesh refinement into polyTopoChange.
        collapser.setRefinement(meshMod);

        // Do all changes
        Info<< "Morphing ..." << endl;

        autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

        collapser.updateMesh(morphMap());

        if (morphMap().hasMotionPoints())
        {
            mesh.movePoints(morphMap().preMotionPoints());
        }

        meshChanged = true;
    }

    if (meshChanged)
    {
        // Write resulting mesh
        if (!overwrite)
        {
            runTime++;
        }
        else
        {
            mesh.setInstance(oldInstance);
        }

        Info<< "Writing collapsed mesh to time " << runTime.timeName() << endl;

        mesh.write();
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
