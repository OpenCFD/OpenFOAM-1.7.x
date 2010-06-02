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
    Extracts triSurface from a polyMesh. Triangulates all boundary faces.
    Region numbers on triangles are the patch numbers of the polyMesh.
    Optionally only triangulates named patches.

    If run in parallel the processor patches get filtered out by default and
    the mesh gets merged. (based on vertex position, not topology, so might go
    wrong!).

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "processorPolyPatch.H"
#include "ListListOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("output file");
#   include "addRegionOption.H"
    argList::validOptions.insert("excludeProcPatches", "");
    argList::validOptions.insert("patches", "(patch0 .. patchN)");

#   include "setRootCase.H"
#   include "createTime.H"

    fileName outFileName(runTime.path()/args.additionalArgs()[0]);

    Info<< "Extracting triSurface from boundaryMesh ..."
        << endl << endl;

    Pout<< "Reading mesh from time " << runTime.value() << endl;

#   include "createNamedPolyMesh.H"

    bool includeProcPatches =
       !(
            args.optionFound("excludeProcPatches")
         || Pstream::parRun()
        );

    // Create local surface from:
    // - explicitly named patches only (-patches option)
    // - all patches (default in sequential mode)
    // - all non-processor patches (default in parallel mode)
    // - all non-processor patches (sequential mode, -excludeProcPatches option)

    // Construct table of patches to include.
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    labelHashSet includePatches(bMesh.size());

    if (args.optionFound("patches"))
    {
        wordList patchNames(args.optionLookup("patches")());

        forAll(patchNames, patchNameI)
        {
            const word& patchName = patchNames[patchNameI];

            label patchI = bMesh.findPatchID(patchName);

            if (patchI == -1)
            {
                FatalErrorIn(args.executable()) << "No such patch "
                    << patchName << endl << "Patches are " << bMesh.names()
                    << exit(FatalError);
            }
            includePatches.insert(patchI);
        }
    }
    else
    {
        forAll(bMesh, patchI)
        {
            const polyPatch& patch = bMesh[patchI];

            if (includeProcPatches || !isA<processorPolyPatch>(patch))
            {
                includePatches.insert(patchI);
            }
            else
            {
                Pout<< patch.name() << " : skipped since processorPatch"
                    << endl;
            }
        }
    }

    triSurface localSurface
    (
        triSurfaceTools::triangulate
        (
            mesh.boundaryMesh(),
            includePatches
        )
    );



    if (!Pstream::parRun())
    {
        Info<< "Writing surface to " << outFileName << endl;

        localSurface.write(outFileName);
    }
    else
    {
        // Write local surface
        fileName localPath = runTime.path()/runTime.caseName() + ".ftr";

        Pout<< "Writing local surface to " << localPath << endl;

        localSurface.write(localPath);

        Info<< endl;


        // Gather all points on master
        List<pointField> gatheredPoints(Pstream::nProcs());

        gatheredPoints[Pstream::myProcNo()] = localSurface.points();

        Pstream::gatherList(gatheredPoints);


        // Gather all localSurface patches
        List<geometricSurfacePatchList> gatheredPatches(Pstream::nProcs());

        gatheredPatches[Pstream::myProcNo()] = localSurface.patches();

        Pstream::gatherList(gatheredPatches);


        // Gather all faces
        List<List<labelledTri> > gatheredFaces(Pstream::nProcs());

        gatheredFaces[Pstream::myProcNo()] = localSurface;

        Pstream::gatherList(gatheredFaces);


        if (Pstream::master())
        {
            // On master combine all points
            pointField allPoints =
                ListListOps::combine<pointField>
                (
                    gatheredPoints,
                    accessOp<pointField>()
                );

            // Count number of patches.
            label nPatches = 0;

            forAll(gatheredPatches, procI)
            {
                nPatches += gatheredPatches[procI].size();
            }

            // Count number of faces.
            label nFaces = 0;

            forAll(gatheredFaces, procI)
            {
                nFaces += gatheredFaces[procI].size();
            }



            // Loop over all processors and
            // - construct mapping from local to global patches
            // - relabel faces (both points and regions)

            label newPatchI = 0;

            // Name to new patchI
            HashTable<label> nameToIndex(2*nPatches);

            // Storage (oversized) for all patches
            geometricSurfacePatchList allPatches(nPatches);

            label newFaceI = 0;

            // Storage for all faces
            List<labelledTri> allFaces(nFaces);

            // Offset into allPoints for current processor
            label pointOffset = 0;

            for (label procI = 0; procI < Pstream::nProcs(); procI++)
            {
                Info<< "Processor " << procI << endl
                    << "-----------" << endl;

                const geometricSurfacePatchList& patches =
                    gatheredPatches[procI];

                // From local patch numbering to global
                labelList localToGlobal(patches.size());

                forAll(patches, patchI)
                {
                    const geometricSurfacePatch& sp = patches[patchI];

                    if (!nameToIndex.found(sp.name()))
                    {
                        nameToIndex.insert(sp.name(), newPatchI);

                        localToGlobal[patchI] = newPatchI;

                        allPatches[newPatchI] = sp;

                        newPatchI++;
                    }
                    else
                    {
                        localToGlobal[patchI] = nameToIndex[sp.name()];
                    }
                }

                Info<< "Local patch to global patch mapping:"
                    << endl;

                forAll(patches, patchI)
                {
                    Info<< "    name   : " << patches[patchI].name() << endl
                        << "    local  : " << patchI << endl
                        << "    global : " << localToGlobal[patchI]
                        << endl;
                }

                Info<< "Local points added in global points starting at "
                    << pointOffset
                    << endl;


                // Collect and relabel faces
                const List<labelledTri>& localFaces = gatheredFaces[procI];


                forAll(localFaces, faceI)
                {
                    const labelledTri& f = localFaces[faceI];

                    allFaces[newFaceI++] =
                        labelledTri
                        (
                            f[0] + pointOffset,
                            f[1] + pointOffset,
                            f[2] + pointOffset,
                            localToGlobal[f.region()]
                        );
                }

                pointOffset += gatheredPoints[procI].size();

                Info<< endl;
            }
            allPatches.setSize(newPatchI);

            // We now have allPoints, allFaces and allPatches.
            // Construct overall (yet unmerged) surface from these.

            triSurface allSurf(allFaces, allPatches, allPoints);

            // Cleanup (which does point merge as well
            allSurf.cleanup(false);

            // Write surface mesh

            label slashIndex = runTime.caseName().find_last_of('/');

            fileName globalCasePath(runTime.caseName().substr(0, slashIndex));

            Info<< "Writing merged surface to " << globalCasePath << endl;

            // create database for the sequential run
            fileName globalPath
            (
                runTime.rootPath()
              / globalCasePath
              / args.additionalArgs()[0]
            );

            allSurf.write(globalPath);
        }
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
