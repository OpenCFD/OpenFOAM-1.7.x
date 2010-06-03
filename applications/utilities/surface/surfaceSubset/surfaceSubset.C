/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    A surface analysis tool which sub-sets the triSurface
    to choose only a part of interest. Based on subsetMesh.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "argList.H"
#include "OFstream.H"
#include "IFstream.H"
#include "Switch.H"
#include "IOdictionary.H"
#include "boundBox.H"
#include "indexedOctree.H"
#include "octree.H"
#include "treeDataTriSurface.H"
#include "Random.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("surfaceSubsetDict");
    argList::validArgs.append("surface file");
    argList::validArgs.append("output file");
    argList args(argc, argv);

    Info<< "Reading dictionary " << args.additionalArgs()[0] << " ..." << endl;
    IFstream dictFile(args.additionalArgs()[0]);
    dictionary meshSubsetDict(dictFile);

    Info<< "Reading surface " << args.additionalArgs()[1] << " ..." << endl;
    triSurface surf1(args.additionalArgs()[1]);

    Info<< "Original:" << endl;
    surf1.writeStats(Info);
    Info<< endl;


    labelList markedPoints
    (
        meshSubsetDict.lookup("localPoints")
    );

    labelList markedEdges
    (
        meshSubsetDict.lookup("edges")
    );

    labelList markedFaces
    (
        meshSubsetDict.lookup("faces")
    );

    pointField markedZone
    (
        meshSubsetDict.lookup("zone")
    );

    if (markedZone.size() && markedZone.size() != 2)
    {
        FatalErrorIn(args.executable())
            << "zone specification should be two points, min and max of "
            << "the boundingbox" << endl
            << "zone:" << markedZone
            << exit(FatalError);
    }

    Switch addFaceNeighbours
    (
        meshSubsetDict.lookup("addFaceNeighbours")
    );

    Switch invertSelection
    (
        meshSubsetDict.lookup("invertSelection")
    );

    // Mark the cells for the subset

    // Faces to subset
    boolList facesToSubset(surf1.size(), false);


    //
    // pick up faces connected to "localPoints"
    //

    if (markedPoints.size())
    {
        Info << "Found " << markedPoints.size() << " marked point(s)." << endl;

        // pick up cells sharing the point

        forAll (markedPoints, pointI)
        {
            if
            (
                markedPoints[pointI] < 0
             || markedPoints[pointI] >= surf1.nPoints()
            )
            {
                FatalErrorIn(args.executable())
                    << "localPoint label " << markedPoints[pointI]
                    << "out of range."
                    << " The mesh has got "
                    << surf1.nPoints() << " localPoints."
                    << exit(FatalError);
            }

            const labelList& curFaces =
                surf1.pointFaces()[markedPoints[pointI]];

            forAll (curFaces, i)
            {
                facesToSubset[curFaces[i]] =  true;
            }
        }
    }



    //
    // pick up faces connected to "edges"
    //

    if (markedEdges.size())
    {
        Info << "Found " << markedEdges.size() << " marked edge(s)." << endl;

        // pick up cells sharing the edge

        forAll (markedEdges, edgeI)
        {
            if
            (
                markedEdges[edgeI] < 0
             || markedEdges[edgeI] >= surf1.nEdges()
            )
            {
                FatalErrorIn(args.executable())
                    << "edge label " << markedEdges[edgeI]
                    << "out of range."
                    << " The mesh has got "
                    << surf1.nEdges() << " edges."
                    << exit(FatalError);
            }

            const labelList& curFaces = surf1.edgeFaces()[markedEdges[edgeI]];

            forAll (curFaces, i)
            {
                facesToSubset[curFaces[i]] =  true;
            }
        }
    }


    //
    // pick up faces with centre inside "zone"
    //

    if (markedZone.size() == 2)
    {
        const point& min = markedZone[0];
        const point& max = markedZone[1];

        Info << "Using zone min:" << min << " max:" << max << endl;

        forAll(surf1, faceI)
        {
            const labelledTri& f = surf1[faceI];
            const point centre = f.centre(surf1.points());

            if
            (
                (centre.x() >= min.x())
             && (centre.y() >= min.y())
             && (centre.z() >= min.z())
             && (centre.x() <= max.x())
             && (centre.y() <= max.y())
             && (centre.z() <= max.z())
            )
            {
                facesToSubset[faceI] = true;
            }
        }
    }


    //
    // pick up faces on certain side of surface
    //

    if (meshSubsetDict.found("surface"))
    {
        const dictionary& surfDict = meshSubsetDict.subDict("surface");

        fileName surfName(surfDict.lookup("name"));

        Switch outside(surfDict.lookup("outside"));

        if (outside)
        {
            Info<< "Selecting all triangles with centre outside surface "
                << surfName << endl;
        }
        else
        {
            Info<< "Selecting all triangles with centre inside surface "
                << surfName << endl;
        }

        // Read surface to select on
        triSurface selectSurf(surfName);

        // bb of surface
        treeBoundBox bb(selectSurf.localPoints());

        // Radnom number generator
        Random rndGen(354543);

        // search engine
        indexedOctree<treeDataTriSurface> selectTree
        (
            treeDataTriSurface(selectSurf),
            bb.extend(rndGen, 1E-4),    // slightly randomize bb
            8,      // maxLevel
            10,     // leafsize
            3.0     // duplicity
        );

        // Check if face (centre) is in outside or inside.
        forAll(facesToSubset, faceI)
        {
            if (!facesToSubset[faceI])
            {
                const point fc(surf1[faceI].centre(surf1.points()));

                indexedOctree<treeDataTriSurface>::volumeType t =
                    selectTree.getVolumeType(fc);

                if (t == indexedOctree<treeDataTriSurface>::INSIDE && !outside)
                {
                    facesToSubset[faceI] = true;
                }
                else if
                (
                    t == indexedOctree<treeDataTriSurface>::OUTSIDE
                 && outside
                )
                {
                    facesToSubset[faceI] = true;
                }
            }
        }
    }


    //
    // pick up specified "faces"
    //

    // Number of additional faces picked up because of addFaceNeighbours
    label nFaceNeighbours = 0;

    if (markedFaces.size())
    {
        Info << "Found " << markedFaces.size() << " marked face(s)." << endl;

        // Check and mark faces to pick up
        forAll (markedFaces, faceI)
        {
            if
            (
                markedFaces[faceI] < 0
             || markedFaces[faceI] >= surf1.size()
            )
            {
                FatalErrorIn(args.executable())
                    << "Face label " << markedFaces[faceI] << "out of range."
                    << " The mesh has got "
                    << surf1.size() << " faces."
                    << exit(FatalError);
            }

            // Mark the face
            facesToSubset[markedFaces[faceI]] = true;

            // mark its neighbours if requested
            if (addFaceNeighbours)
            {
                const labelList& curFaces =
                    surf1.faceFaces()[markedFaces[faceI]];

                forAll (curFaces, i)
                {
                    label faceI = curFaces[i];

                    if (!facesToSubset[faceI])
                    {
                        facesToSubset[faceI] =  true;
                        nFaceNeighbours++;
                    }
                }
            }
        }
    }

    if (addFaceNeighbours)
    {
        Info<< "Added " << nFaceNeighbours
            << " faces because of addFaceNeighbours" << endl;
    }


    if (invertSelection)
    {
        Info<< "Inverting selection." << endl;
        boolList newFacesToSubset(facesToSubset.size());

        forAll(facesToSubset, i)
        {
            if (facesToSubset[i])
            {
                newFacesToSubset[i] = false;
            }
            else
            {
                newFacesToSubset[i] = true;
            }
        }
        facesToSubset.transfer(newFacesToSubset);
    }


    // Create subsetted surface
    labelList pointMap;
    labelList faceMap;
    triSurface surf2
    (
        surf1.subsetMesh(facesToSubset, pointMap, faceMap)
    );

    Info<< "Subset:" << endl;
    surf2.writeStats(Info);
    Info << endl;

    fileName outFileName(args.additionalArgs()[2]);

    Info << "Writing surface to " << outFileName << endl;

    surf2.write(outFileName);

    return 0;
}


// ************************************************************************* //
