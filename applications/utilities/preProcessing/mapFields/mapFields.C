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

Application
    mapFields

Description
    Maps volume fields from one mesh to another, reading and
    interpolating all fields present in the time directory of both cases.
    Parallel and non-parallel cases are handled without the need to reconstruct
    them first.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "meshToMesh.H"
#include "MapVolFields.H"
#include "MapConsistentVolFields.H"
#include "UnMapped.H"
#include "processorFvPatch.H"
#include "mapLagrangian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mapConsistentMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget
)
{
    // Create the interpolation scheme
    meshToMesh meshToMeshInterp(meshSource, meshTarget);

    Info<< nl
        << "Consistently creating and mapping fields for time "
        << meshSource.time().timeName() << nl << endl;

    {
        // Search for list of objects for this time
        IOobjectList objects(meshSource, meshSource.time().timeName());

        // Map volFields
        // ~~~~~~~~~~~~~
        MapConsistentVolFields<scalar>(objects, meshToMeshInterp);
        MapConsistentVolFields<vector>(objects, meshToMeshInterp);
        MapConsistentVolFields<sphericalTensor>(objects, meshToMeshInterp);
        MapConsistentVolFields<symmTensor>(objects, meshToMeshInterp);
        MapConsistentVolFields<tensor>(objects, meshToMeshInterp);
    }

    {
        // Search for list of target objects for this time
        IOobjectList objects(meshTarget, meshTarget.time().timeName());

        // Mark surfaceFields as unmapped
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        UnMapped<surfaceScalarField>(objects);
        UnMapped<surfaceVectorField>(objects);
        UnMapped<surfaceSphericalTensorField>(objects);
        UnMapped<surfaceSymmTensorField>(objects);
        UnMapped<surfaceTensorField>(objects);

        // Mark pointFields as unmapped
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        UnMapped<pointScalarField>(objects);
        UnMapped<pointVectorField>(objects);
        UnMapped<pointSphericalTensorField>(objects);
        UnMapped<pointSymmTensorField>(objects);
        UnMapped<pointTensorField>(objects);
    }

    mapLagrangian(meshToMeshInterp);
}


void mapSubMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatches
)
{
    // Create the interpolation scheme
    meshToMesh meshToMeshInterp
    (
        meshSource,
        meshTarget,
        patchMap,
        cuttingPatches
    );

    Info<< nl
        << "Mapping fields for time " << meshSource.time().timeName()
        << nl << endl;

    {
        // Search for list of source objects for this time
        IOobjectList objects(meshSource, meshSource.time().timeName());

        // Map volFields
        // ~~~~~~~~~~~~~
        MapVolFields<scalar>(objects, meshToMeshInterp);
        MapVolFields<vector>(objects, meshToMeshInterp);
        MapVolFields<sphericalTensor>(objects, meshToMeshInterp);
        MapVolFields<symmTensor>(objects, meshToMeshInterp);
        MapVolFields<tensor>(objects, meshToMeshInterp);
    }

    {
        // Search for list of target objects for this time
        IOobjectList objects(meshTarget, meshTarget.time().timeName());

        // Mark surfaceFields as unmapped
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        UnMapped<surfaceScalarField>(objects);
        UnMapped<surfaceVectorField>(objects);
        UnMapped<surfaceSphericalTensorField>(objects);
        UnMapped<surfaceSymmTensorField>(objects);
        UnMapped<surfaceTensorField>(objects);

        // Mark pointFields as unmapped
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        UnMapped<pointScalarField>(objects);
        UnMapped<pointVectorField>(objects);
        UnMapped<pointSphericalTensorField>(objects);
        UnMapped<pointSymmTensorField>(objects);
        UnMapped<pointTensorField>(objects);
    }

    mapLagrangian(meshToMeshInterp);
}


void mapConsistentSubMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget
)
{
    HashTable<word> patchMap;
    HashTable<label> cuttingPatchTable;

    forAll(meshTarget.boundary(), patchi)
    {
        if (typeid(meshTarget.boundary()[patchi]) != typeid(processorFvPatch))
        {
            patchMap.insert
            (
                meshTarget.boundary()[patchi].name(),
                meshTarget.boundary()[patchi].name()
            );
        }
        else
        {
            cuttingPatchTable.insert
            (
                meshTarget.boundaryMesh()[patchi].name(),
                -1
            );
        }
    }

    mapSubMesh(meshSource, meshTarget, patchMap, cuttingPatchTable.toc());
}


wordList addProcessorPatches
(
    const fvMesh& meshTarget,
    const wordList& cuttingPatches
)
{
    // Add the processor patches to the cutting list
    HashTable<label> cuttingPatchTable;
    forAll (cuttingPatches, i)
    {
        cuttingPatchTable.insert(cuttingPatches[i], i);
    }

    forAll (meshTarget.boundary(), patchi)
    {
        if (typeid(meshTarget.boundary()[patchi]) == typeid(processorFvPatch))
        {
            if
            (
               !cuttingPatchTable.found
                (
                    meshTarget.boundaryMesh()[patchi].name()
                )
            )
            {
                cuttingPatchTable.insert
                (
                    meshTarget.boundaryMesh()[patchi].name(),
                    -1
                );
            }
        }
    }

    return cuttingPatchTable.toc();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRoots.H"
    #include "createTimes.H"

    HashTable<word> patchMap;
    wordList cuttingPatches;

    if (!consistent)
    {
        IOdictionary mapFieldsDict
        (
            IOobject
            (
                "mapFieldsDict",
                runTimeTarget.system(),
                runTimeTarget,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        mapFieldsDict.lookup("patchMap") >> patchMap;

        mapFieldsDict.lookup("cuttingPatches") >>  cuttingPatches;
    }

    if (parallelSource && !parallelTarget)
    {
        IOdictionary decompositionDict
        (
            IOobject
            (
                "decomposeParDict",
                runTimeSource.system(),
                runTimeSource,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        int nProcs(readInt(decompositionDict.lookup("numberOfSubdomains")));

        Info<< "Create target mesh\n" << endl;

        fvMesh meshTarget
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTimeTarget.timeName(),
                runTimeTarget
            )
        );

        Info<< "Target mesh size: " << meshTarget.nCells() << endl;

        for (int procI=0; procI<nProcs; procI++)
        {
            Info<< nl << "Source processor " << procI << endl;

            Time runTimeSource
            (
                Time::controlDictName,
                rootDirSource,
                caseDirSource/fileName(word("processor") + name(procI))
            );

            #include "setTimeIndex.H"

            fvMesh meshSource
            (
                IOobject
                (
                    fvMesh::defaultRegion,
                    runTimeSource.timeName(),
                    runTimeSource
                )
            );

            Info<< "mesh size: " << meshSource.nCells() << endl;

            if (consistent)
            {
                mapConsistentSubMesh(meshSource, meshTarget);
            }
            else
            {
                mapSubMesh(meshSource, meshTarget, patchMap, cuttingPatches);
            }
        }
    }
    else if (!parallelSource && parallelTarget)
    {
        IOdictionary decompositionDict
        (
            IOobject
            (
                "decomposeParDict",
                runTimeTarget.system(),
                runTimeTarget,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        int nProcs(readInt(decompositionDict.lookup("numberOfSubdomains")));

        Info<< "Create source mesh\n" << endl;

        #include "setTimeIndex.H"

        fvMesh meshSource
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTimeSource.timeName(),
                runTimeSource
            )
        );

        Info<< "Source mesh size: " << meshSource.nCells() << endl;

        for (int procI=0; procI<nProcs; procI++)
        {
            Info<< nl << "Target processor " << procI << endl;

            Time runTimeTarget
            (
                Time::controlDictName,
                rootDirTarget,
                caseDirTarget/fileName(word("processor") + name(procI))
            );

            fvMesh meshTarget
            (
                IOobject
                (
                    fvMesh::defaultRegion,
                    runTimeTarget.timeName(),
                    runTimeTarget
                )
            );

            Info<< "mesh size: " << meshTarget.nCells() << endl;

            if (consistent)
            {
                mapConsistentSubMesh(meshSource, meshTarget);
            }
            else
            {
                mapSubMesh
                (
                    meshSource,
                    meshTarget,
                    patchMap,
                    addProcessorPatches(meshTarget, cuttingPatches)
                );
            }
        }
    }
    else if (parallelSource && parallelTarget)
    {
        IOdictionary decompositionDictSource
        (
            IOobject
            (
                "decomposeParDict",
                runTimeSource.system(),
                runTimeSource,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        int nProcsSource
        (
            readInt(decompositionDictSource.lookup("numberOfSubdomains"))
        );


        IOdictionary decompositionDictTarget
        (
            IOobject
            (
                "decomposeParDict",
                runTimeTarget.system(),
                runTimeTarget,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        int nProcsTarget
        (
            readInt(decompositionDictTarget.lookup("numberOfSubdomains"))
        );

        List<boundBox> bbsTarget(nProcsTarget);
        List<bool> bbsTargetSet(nProcsTarget, false);

        for (int procISource=0; procISource<nProcsSource; procISource++)
        {
            Info<< nl << "Source processor " << procISource << endl;

            Time runTimeSource
            (
                Time::controlDictName,
                rootDirSource,
                caseDirSource/fileName(word("processor") + name(procISource))
            );

            #include "setTimeIndex.H"

            fvMesh meshSource
            (
                IOobject
                (
                    fvMesh::defaultRegion,
                    runTimeSource.timeName(),
                    runTimeSource
                )
            );

            Info<< "mesh size: " << meshSource.nCells() << endl;

            boundBox bbSource(meshSource.bounds());

            for (int procITarget=0; procITarget<nProcsTarget; procITarget++)
            {
                if
                (
                    !bbsTargetSet[procITarget]
                  || (
                      bbsTargetSet[procITarget]
                   && bbsTarget[procITarget].overlaps(bbSource)
                     )
                )
                {
                    Info<< nl << "Target processor " << procITarget << endl;

                    Time runTimeTarget
                    (
                        Time::controlDictName,
                        rootDirTarget,
                        caseDirTarget/fileName(word("processor")
                      + name(procITarget))
                    );

                    fvMesh meshTarget
                    (
                        IOobject
                        (
                            fvMesh::defaultRegion,
                            runTimeTarget.timeName(),
                            runTimeTarget
                        )
                    );

                    Info<< "mesh size: " << meshTarget.nCells() << endl;

                    bbsTarget[procITarget] = meshTarget.bounds();
                    bbsTargetSet[procITarget] = true;

                    if (bbsTarget[procITarget].overlaps(bbSource))
                    {
                        if (consistent)
                        {
                            mapConsistentSubMesh(meshSource, meshTarget);
                        }
                        else
                        {
                            mapSubMesh
                            (
                                meshSource,
                                meshTarget,
                                patchMap,
                                addProcessorPatches(meshTarget, cuttingPatches)
                            );
                        }
                    }
                }
            }
        }
    }
    else
    {
        #include "setTimeIndex.H"

        Info<< "Create meshes\n" << endl;

        fvMesh meshSource
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTimeSource.timeName(),
                runTimeSource
            )
        );

        fvMesh meshTarget
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTimeTarget.timeName(),
                runTimeTarget
            )
        );

        Info<< "Source mesh size: " << meshSource.nCells() << tab
            << "Target mesh size: " << meshTarget.nCells() << nl << endl;

        if (consistent)
        {
            mapConsistentMesh(meshSource, meshTarget);
        }
        else
        {
            mapSubMesh(meshSource, meshTarget, patchMap, cuttingPatches);
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
