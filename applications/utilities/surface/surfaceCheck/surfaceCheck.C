/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "triangle.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "argList.H"
#include "OFstream.H"
#include "surfaceIntersection.H"
#include "SortableList.H"
#include "PatchTools.H"

using namespace Foam;

// Does face use valid vertices?
bool validTri(const bool verbose, const triSurface& surf, const label faceI)
{
    // Simple check on indices ok.

    const labelledTri& f = surf[faceI];

    if
    (
        (f[0] < 0) || (f[0] >= surf.points().size())
     || (f[1] < 0) || (f[1] >= surf.points().size())
     || (f[2] < 0) || (f[2] >= surf.points().size())
    )
    {
        WarningIn("validTri(const triSurface&, const label)")
            << "triangle " << faceI << " vertices " << f
            << " uses point indices outside point range 0.."
            << surf.points().size()-1 << endl;

        return false;
    }

    if ((f[0] == f[1]) || (f[0] == f[2]) || (f[1] == f[2]))
    {
        WarningIn("validTri(const triSurface&, const label)")
            << "triangle " << faceI
            << " uses non-unique vertices " << f
            << " coords:" << f.points(surf.points())
            << endl;
        return false;
    }

    // duplicate triangle check

    const labelList& fFaces = surf.faceFaces()[faceI];

    // Check if faceNeighbours use same points as this face.
    // Note: discards normal information - sides of baffle are merged.
    forAll(fFaces, i)
    {
        label nbrFaceI = fFaces[i];

        if (nbrFaceI <= faceI)
        {
            // lower numbered faces already checked
            continue;
        }

        const labelledTri& nbrF = surf[nbrFaceI];

        if
        (
            ((f[0] == nbrF[0]) || (f[0] == nbrF[1]) || (f[0] == nbrF[2]))
         && ((f[1] == nbrF[0]) || (f[1] == nbrF[1]) || (f[1] == nbrF[2]))
         && ((f[2] == nbrF[0]) || (f[2] == nbrF[1]) || (f[2] == nbrF[2]))
        )
        {
            WarningIn("validTri(const triSurface&, const label)")
                << "triangle " << faceI << " vertices " << f
                << " has the same vertices as triangle " << nbrFaceI
                << " vertices " << nbrF
                << " coords:" << f.points(surf.points())
                << endl;

            return false;
        }
    }
    return true;
}


labelList countBins
(
    const scalar min,
    const scalar max,
    const label nBins,
    const scalarField& vals
)
{
    scalar dist = nBins/(max - min);

    labelList binCount(nBins, 0);

    forAll(vals, i)
    {
        scalar val = vals[i];

        label index = -1;

        if (Foam::mag(val - min) < SMALL)
        {
            index = 0;
        }
        else if (val >= max - SMALL)
        {
            index = nBins - 1;
        }
        else
        {
            index = label((val - min)*dist);

            if ((index < 0) || (index >= nBins))
            {
                WarningIn
                (
                    "countBins(const scalar, const scalar, const label"
                    ", const scalarField&)"
                )   << "value " << val << " at index " << i
                    << " outside range " << min << " .. " << max << endl;

                if (index < 0)
                {
                    index = 0;
                }
                else
                {
                    index = nBins - 1;
                }
            }
        }
        binCount[index]++;
    }

    return binCount;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();

    argList::validArgs.clear();
    argList::validArgs.append("surface file");
    argList::validOptions.insert("checkSelfIntersection", "");
    argList::validOptions.insert("verbose", "");
    argList args(argc, argv);

    bool checkSelfIntersection = args.optionFound("checkSelfIntersection");
    bool verbose = args.optionFound("verbose");

    fileName surfFileName(args.additionalArgs()[0]);
    Pout<< "Reading surface from " << surfFileName << " ..." << nl << endl;


    // Read
    // ~~~~

    triSurface surf(surfFileName);


    Pout<< "Statistics:" << endl;
    surf.writeStats(Pout);
    Pout<< endl;


    // Region sizes
    // ~~~~~~~~~~~~

    {
        labelList regionSize(surf.patches().size(), 0);

        forAll(surf, faceI)
        {
            label region = surf[faceI].region();

            if (region < 0 || region >= regionSize.size())
            {
                WarningIn(args.executable())
                    << "Triangle " << faceI << " vertices " << surf[faceI]
                    << " has region " << region << " which is outside the range"
                    << " of regions 0.." << surf.patches().size()-1
                    << endl;
            }
            else
            {
                regionSize[region]++;
            }
        }

        Pout<< "Region\tSize" << nl
            << "------\t----" << nl;
        forAll(surf.patches(), patchI)
        {
            Pout<< surf.patches()[patchI].name() << '\t'
                << regionSize[patchI] << nl;
        }
        Pout<< nl << endl;
    }


    // Check triangles
    // ~~~~~~~~~~~~~~~

    {
        DynamicList<label> illegalFaces(surf.size()/100 + 1);

        forAll(surf, faceI)
        {
            if (!validTri(verbose, surf, faceI))
            {
                illegalFaces.append(faceI);
            }
        }

        if (illegalFaces.size())
        {
            Pout<< "Surface has " << illegalFaces.size()
                << " illegal triangles." << endl;

            OFstream str("illegalFaces");
            Pout<< "Dumping conflicting face labels to " << str.name() << endl
                << "Paste this into the input for surfaceSubset" << endl;
            str << illegalFaces;
        }
        else
        {
            Pout<< "Surface has no illegal triangles." << endl;
        }
        Pout<< endl;
    }



    // Triangle quality
    // ~~~~~~~~~~~~~~~~

    {
        scalarField triQ(surf.size(), 0);
        forAll(surf, faceI)
        {
            const labelledTri& f = surf[faceI];

            if (f[0] == f[1] || f[0] == f[2] || f[1] == f[2])
            {
                //WarningIn(args.executable())
                //    << "Illegal triangle " << faceI << " vertices " << f
                //    << " coords " << f.points(surf.points()) << endl;
            }
            else
            {
                triPointRef tri
                (
                    surf.points()[f[0]],
                    surf.points()[f[1]],
                    surf.points()[f[2]]
                );

                vector ba(tri.b() - tri.a());
                ba /= mag(ba) + VSMALL;

                vector ca(tri.c() - tri.a());
                ca /= mag(ca) + VSMALL;

                if (mag(ba&ca) > 1-1E-3)
                {
                    triQ[faceI] = SMALL;
                }
                else
                {
                    triQ[faceI] = triPointRef
                    (
                        surf.points()[f[0]],
                        surf.points()[f[1]],
                        surf.points()[f[2]]
                    ).quality();
                }
            }
        }

        labelList binCount = countBins(0, 1, 20, triQ);

        Pout<< "Triangle quality (equilateral=1, collapsed=0):"
            << endl;


        OSstream& os = Pout;
        os.width(4);

        scalar dist = (1.0 - 0.0)/20.0;
        scalar min = 0;
        forAll(binCount, binI)
        {
            Pout<< "    " << min << " .. " << min+dist << "  : "
                << 1.0/surf.size() * binCount[binI]
                << endl;
            min += dist;
        }
        Pout<< endl;

        label minIndex = findMin(triQ);
        label maxIndex = findMax(triQ);

        Pout<< "    min " << triQ[minIndex] << " for triangle " << minIndex
            << nl
            << "    max " << triQ[maxIndex] << " for triangle " << maxIndex
            << nl
            << endl;


        if (triQ[minIndex] < SMALL)
        {
            WarningIn(args.executable()) << "Minimum triangle quality is "
                << triQ[minIndex] << ". This might give problems in"
                << " self-intersection testing later on." << endl;
        }

        // Dump for subsetting
        {
            DynamicList<label> problemFaces(surf.size()/100+1);

            forAll(triQ, faceI)
            {
                if (triQ[faceI] < 1E-11)
                {
                    problemFaces.append(faceI);
                }
            }
            OFstream str("badFaces");

            Pout<< "Dumping bad quality faces to " << str.name() << endl
                << "Paste this into the input for surfaceSubset" << nl
                << nl << endl;

            str << problemFaces;
        }
    }



    // Edges
    // ~~~~~
    {
        const edgeList& edges = surf.edges();
        const pointField& localPoints = surf.localPoints();

        scalarField edgeMag(edges.size());

        forAll(edges, edgeI)
        {
            edgeMag[edgeI] = edges[edgeI].mag(localPoints);
        }

        label minEdgeI = findMin(edgeMag);
        label maxEdgeI = findMax(edgeMag);

        const edge& minE = edges[minEdgeI];
        const edge& maxE = edges[maxEdgeI];


        Pout<< "Edges:" << nl
            << "    min " << edgeMag[minEdgeI] << " for edge " << minEdgeI
            << " points " << localPoints[minE[0]] << localPoints[minE[1]]
            << nl
            << "    max " << edgeMag[maxEdgeI] << " for edge " << maxEdgeI
            << " points " << localPoints[maxE[0]] << localPoints[maxE[1]]
            << nl
            << endl;
    }



    // Close points
    // ~~~~~~~~~~~~
    {
        const edgeList& edges = surf.edges();
        const pointField& localPoints = surf.localPoints();

        const boundBox bb(localPoints);
        scalar smallDim = 1E-6 * bb.mag();

        Pout<< "Checking for points less than 1E-6 of bounding box ("
            << bb.span() << " meter) apart."
            << endl;

        // Sort points
        SortableList<scalar> sortedMag(mag(localPoints));

        label nClose = 0;

        for (label i = 1; i < sortedMag.size(); i++)
        {
            label ptI = sortedMag.indices()[i];

            label prevPtI = sortedMag.indices()[i-1];

            if (mag(localPoints[ptI] - localPoints[prevPtI]) < smallDim)
            {
                // Check if neighbours.
                const labelList& pEdges = surf.pointEdges()[ptI];

                label edgeI = -1;

                forAll(pEdges, i)
                {
                    const edge& e = edges[pEdges[i]];

                    if (e[0] == prevPtI || e[1] == prevPtI)
                    {
                        // point1 and point0 are connected through edge.
                        edgeI = pEdges[i];

                        break;
                    }
                }

                nClose++;

                if (edgeI == -1)
                {
                    Pout<< "    close unconnected points "
                        << ptI << ' ' << localPoints[ptI]
                        << " and " << prevPtI << ' '
                        << localPoints[prevPtI]
                        << " distance:"
                        << mag(localPoints[ptI] - localPoints[prevPtI])
                        << endl;
                }
                else
                {
                    Pout<< "    small edge between points "
                        << ptI << ' ' << localPoints[ptI]
                        << " and " << prevPtI << ' '
                        << localPoints[prevPtI]
                        << " distance:"
                        << mag(localPoints[ptI] - localPoints[prevPtI])
                        << endl;
                }
            }
        }

        Pout<< "Found " << nClose << " nearby points." << nl
            << endl;
    }



    // Check manifold
    // ~~~~~~~~~~~~~~

    DynamicList<label> problemFaces(surf.size()/100 + 1);

    const labelListList& eFaces = surf.edgeFaces();

    label nSingleEdges = 0;
    forAll(eFaces, edgeI)
    {
        const labelList& myFaces = eFaces[edgeI];

        if (myFaces.size() == 1)
        {
            problemFaces.append(myFaces[0]);

            nSingleEdges++;
        }
    }

    label nMultEdges = 0;
    forAll(eFaces, edgeI)
    {
        const labelList& myFaces = eFaces[edgeI];

        if (myFaces.size() > 2)
        {
            forAll(myFaces, myFaceI)
            {
                problemFaces.append(myFaces[myFaceI]);
            }

            nMultEdges++;
        }
    }
    problemFaces.shrink();

    if ((nSingleEdges != 0) || (nMultEdges != 0))
    {
        Pout<< "Surface is not closed since not all edges connected to "
            << "two faces:" << endl
            << "    connected to one face : " << nSingleEdges << endl
            << "    connected to >2 faces : " << nMultEdges << endl;

        Pout<< "Conflicting face labels:" << problemFaces.size() << endl;

        OFstream str("problemFaces");

        Pout<< "Dumping conflicting face labels to " << str.name() << endl
            << "Paste this into the input for surfaceSubset" << endl;

        str << problemFaces;
    }
    else
    {
        Pout<< "Surface is closed. All edges connected to two faces." << endl;
    }
    Pout<< endl;



    // Check singly connected domain
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList faceZone;
    label numZones = surf.markZones(boolList(surf.nEdges(), false), faceZone);

    Pout<< "Number of unconnected parts : " << numZones << endl;

    if (numZones > 1)
    {
        Pout<< "Splitting surface into parts ..." << endl << endl;

        fileName surfFileNameBase(surfFileName.name());

        for(label zone = 0; zone < numZones; zone++)
        {
            boolList includeMap(surf.size(), false);

            forAll(faceZone, faceI)
            {
                if (faceZone[faceI] == zone)
                {
                    includeMap[faceI] = true;
                }
            }

            labelList pointMap;
            labelList faceMap;

            triSurface subSurf
            (
                surf.subsetMesh
                (
                    includeMap,
                    pointMap,
                    faceMap
                )
            );

            fileName subFileName
            (
                surfFileNameBase.lessExt()
              + "_"
              + name(zone)
              + ".ftr"
            );

            Pout<< "writing part " << zone << " size " << subSurf.size()
                << " to " << subFileName << endl;

            subSurf.write(subFileName);
        }

        return 0;
    }



    // Check orientation
    // ~~~~~~~~~~~~~~~~~

    labelHashSet borderEdge(surf.size()/1000);
    PatchTools::checkOrientation(surf, false, &borderEdge);

    //
    // Colour all faces into zones using borderEdge
    //
    labelList normalZone;
    label numNormalZones = PatchTools::markZones(surf, borderEdge, normalZone);

    Pout<< endl
        << "Number of zones (connected area with consistent normal) : "
        << numNormalZones << endl;

    if (numNormalZones > 1)
    {
        Pout<< "More than one normal orientation." << endl;
    }
    Pout<< endl;



    // Check self-intersection
    // ~~~~~~~~~~~~~~~~~~~~~~~

    if (checkSelfIntersection)
    {
        Pout<< "Checking self-intersection." << endl;

        triSurfaceSearch querySurf(surf);
        surfaceIntersection inter(querySurf);

        if (inter.cutEdges().empty() && inter.cutPoints().empty())
        {
            Pout<< "Surface is not self-intersecting" << endl;
        }
        else
        {
            Pout<< "Surface is self-intersecting" << endl;
            Pout<< "Writing edges of intersection to selfInter.obj" << endl;

            OFstream intStream("selfInter.obj");
            forAll(inter.cutPoints(), cutPointI)
            {
                const point& pt = inter.cutPoints()[cutPointI];

                intStream << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z()
                    << endl;
            }
            forAll(inter.cutEdges(), cutEdgeI)
            {
                const edge& e = inter.cutEdges()[cutEdgeI];

                intStream << "l " << e.start()+1 << ' ' << e.end()+1 << endl;
            }
        }
        Pout<< endl;
    }


    Pout<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
