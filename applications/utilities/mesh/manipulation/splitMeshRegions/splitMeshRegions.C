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

Application
    splitMeshRegions

Description
    Splits mesh into multiple regions.

    Each region is defined as a domain whose cells can all be reached by
    cell-face-cell walking without crossing
    - boundary faces
    - additional faces from faceset (-blockedFaces faceSet).
    - any face inbetween differing cellZones (-cellZones)

    Output is:
    - volScalarField with regions as different scalars (-detectOnly) or
    - mesh with multiple regions or
    - mesh with cells put into cellZones (-makeCellZones)

    Note:
    - cellZonesOnly does not do a walk and uses the cellZones only. Use
    this if you don't mind having disconnected domains in a single region.
    This option requires all cells to be in one (and one only) cellZone.

    - cellZonesFileOnly behaves like -cellZonesOnly but reads the cellZones
    from the specified file. This allows one to explicitly specify the region
    distribution and still have multiple cellZones per region.

    - useCellZonesOnly does not do a walk and uses the cellZones only. Use
    this if you don't mind having disconnected domains in a single region.
    This option requires all cells to be in one (and one only) cellZone.


    - Should work in parallel.
    cellZones can differ on either side of processor boundaries in which case
    the faces get moved from processor patch to directMapped patch. Not
    very well tested.

    - If a cell zone gets split into more than one region it can detect
    the largest matching region (-sloppyCellZones). This will accept any
    region that covers more than 50% of the zone. It has to be a subset
    so cannot have any cells in any other zone.

    - writes maps like decomposePar back to original mesh:
        - pointRegionAddressing : for every point in this region the point in
        the original mesh
        - cellRegionAddressing  :   ,,      cell                ,,  cell    ,,
        - faceRegionAddressing  :   ,,      face                ,,  face in
        the original mesh + 'turning index'. For a face in the same orientation
        this is the original facelabel+1, for a turned face this is -facelabel-1
\*---------------------------------------------------------------------------*/

#include "SortableList.H"
#include "argList.H"
#include "regionSplit.H"
#include "fvMeshSubset.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "faceSet.H"
#include "cellSet.H"
#include "polyTopoChange.H"
#include "removeCells.H"
#include "EdgeMap.H"
#include "syncTools.H"
#include "ReadFields.H"
#include "directMappedWallPolyPatch.H"
#include "zeroGradientFvPatchFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void addPatchFields(fvMesh& mesh, const word& patchFieldType)
{
    HashTable<const GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    for
    (
        typename HashTable<const GeoField*>::const_iterator iter = flds.begin();
        iter != flds.end();
        ++iter
    )
    {
        const GeoField& fld = *iter();

        typename GeoField::GeometricBoundaryField& bfld =
            const_cast<typename GeoField::GeometricBoundaryField&>
            (
                fld.boundaryField()
            );

        label sz = bfld.size();
        bfld.setSize(sz+1);
        bfld.set
        (
            sz,
            GeoField::PatchFieldType::New
            (
                patchFieldType,
                mesh.boundary()[sz],
                fld.dimensionedInternalField()
            )
        );
    }
}


// Remove last patch field
template<class GeoField>
void trimPatchFields(fvMesh& mesh, const label nPatches)
{
    HashTable<const GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    for
    (
        typename HashTable<const GeoField*>::const_iterator iter = flds.begin();
        iter != flds.end();
        ++iter
    )
    {
        const GeoField& fld = *iter();

        const_cast<typename GeoField::GeometricBoundaryField&>
        (
            fld.boundaryField()
        ).setSize(nPatches);
    }
}


// Reorder patch field
template<class GeoField>
void reorderPatchFields(fvMesh& mesh, const labelList& oldToNew)
{
    HashTable<const GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    for
    (
        typename HashTable<const GeoField*>::const_iterator iter = flds.begin();
        iter != flds.end();
        ++iter
    )
    {
        const GeoField& fld = *iter();

        typename GeoField::GeometricBoundaryField& bfld =
            const_cast<typename GeoField::GeometricBoundaryField&>
            (
                fld.boundaryField()
            );

        bfld.reorder(oldToNew);
    }
}


// Adds patch if not yet there. Returns patchID.
label addPatch(fvMesh& mesh, const polyPatch& patch)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());

    label patchI = polyPatches.findPatchID(patch.name());
    if (patchI != -1)
    {
        if (polyPatches[patchI].type() == patch.type())
        {
            // Already there
            return patchI;
        }
        else
        {
            FatalErrorIn("addPatch(fvMesh&, const polyPatch*)")
                << "Already have patch " << patch.name()
                << " but of type " << patch.type()
                << exit(FatalError);
        }
    }


    label insertPatchI = polyPatches.size();
    label startFaceI = mesh.nFaces();

    forAll(polyPatches, patchI)
    {
        const polyPatch& pp = polyPatches[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            insertPatchI = patchI;
            startFaceI = pp.start();
            break;
        }
    }


    // Below is all quite a hack. Feel free to change once there is a better
    // mechanism to insert and reorder patches.

    // Clear local fields and e.g. polyMesh parallelInfo.
    mesh.clearOut();

    label sz = polyPatches.size();

    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    // Add polyPatch at the end
    polyPatches.setSize(sz+1);
    polyPatches.set
    (
        sz,
        patch.clone
        (
            polyPatches,
            insertPatchI,   //index
            0,              //size
            startFaceI      //start
        )
    );
    fvPatches.setSize(sz+1);
    fvPatches.set
    (
        sz,
        fvPatch::New
        (
            polyPatches[sz],  // point to newly added polyPatch
            mesh.boundary()
        )
    );

    addPatchFields<volScalarField>
    (
        mesh,
        calculatedFvPatchField<scalar>::typeName
    );
    addPatchFields<volVectorField>
    (
        mesh,
        calculatedFvPatchField<vector>::typeName
    );
    addPatchFields<volSphericalTensorField>
    (
        mesh,
        calculatedFvPatchField<sphericalTensor>::typeName
    );
    addPatchFields<volSymmTensorField>
    (
        mesh,
        calculatedFvPatchField<symmTensor>::typeName
    );
    addPatchFields<volTensorField>
    (
        mesh,
        calculatedFvPatchField<tensor>::typeName
    );

    // Surface fields

    addPatchFields<surfaceScalarField>
    (
        mesh,
        calculatedFvPatchField<scalar>::typeName
    );
    addPatchFields<surfaceVectorField>
    (
        mesh,
        calculatedFvPatchField<vector>::typeName
    );
    addPatchFields<surfaceSphericalTensorField>
    (
        mesh,
        calculatedFvPatchField<sphericalTensor>::typeName
    );
    addPatchFields<surfaceSymmTensorField>
    (
        mesh,
        calculatedFvPatchField<symmTensor>::typeName
    );
    addPatchFields<surfaceTensorField>
    (
        mesh,
        calculatedFvPatchField<tensor>::typeName
    );

    // Create reordering list
    // patches before insert position stay as is
    labelList oldToNew(sz+1);
    for (label i = 0; i < insertPatchI; i++)
    {
        oldToNew[i] = i;
    }
    // patches after insert position move one up
    for (label i = insertPatchI; i < sz; i++)
    {
        oldToNew[i] = i+1;
    }
    // appended patch gets moved to insert position
    oldToNew[sz] = insertPatchI;

    // Shuffle into place
    polyPatches.reorder(oldToNew);
    fvPatches.reorder(oldToNew);

    reorderPatchFields<volScalarField>(mesh, oldToNew);
    reorderPatchFields<volVectorField>(mesh, oldToNew);
    reorderPatchFields<volSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<volSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<volTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceScalarField>(mesh, oldToNew);
    reorderPatchFields<surfaceVectorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceTensorField>(mesh, oldToNew);

    return insertPatchI;
}


// Reorder and delete patches.
void reorderPatches
(
    fvMesh& mesh,
    const labelList& oldToNew,
    const label nNewPatches
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    // Shuffle into place
    polyPatches.reorder(oldToNew);
    fvPatches.reorder(oldToNew);

    reorderPatchFields<volScalarField>(mesh, oldToNew);
    reorderPatchFields<volVectorField>(mesh, oldToNew);
    reorderPatchFields<volSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<volSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<volTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceScalarField>(mesh, oldToNew);
    reorderPatchFields<surfaceVectorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceTensorField>(mesh, oldToNew);

    // Remove last.
    polyPatches.setSize(nNewPatches);
    fvPatches.setSize(nNewPatches);
    trimPatchFields<volScalarField>(mesh, nNewPatches);
    trimPatchFields<volVectorField>(mesh, nNewPatches);
    trimPatchFields<volSphericalTensorField>(mesh, nNewPatches);
    trimPatchFields<volSymmTensorField>(mesh, nNewPatches);
    trimPatchFields<volTensorField>(mesh, nNewPatches);
    trimPatchFields<surfaceScalarField>(mesh, nNewPatches);
    trimPatchFields<surfaceVectorField>(mesh, nNewPatches);
    trimPatchFields<surfaceSphericalTensorField>(mesh, nNewPatches);
    trimPatchFields<surfaceSymmTensorField>(mesh, nNewPatches);
    trimPatchFields<surfaceTensorField>(mesh, nNewPatches);
}


template<class GeoField>
void subsetVolFields
(
    const fvMesh& mesh,
    const fvMesh& subMesh,
    const labelList& cellMap,
    const labelList& faceMap
)
{
    const labelList patchMap(identity(mesh.boundaryMesh().size()));

    HashTable<const GeoField*> fields
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );
    forAllConstIter(typename HashTable<const GeoField*>, fields, iter)
    {
        const GeoField& fld = *iter();

        Info<< "Mapping field " << fld.name() << endl;

        tmp<GeoField> tSubFld
        (
            fvMeshSubset::interpolate
            (
                fld,
                subMesh,
                patchMap,
                cellMap,
                faceMap
            )
        );

        // Hack: set value to 0 for introduced patches (since don't
        //       get initialised.
        forAll(tSubFld().boundaryField(), patchI)
        {
            const fvPatchField<typename GeoField::value_type>& pfld =
                tSubFld().boundaryField()[patchI];

            if
            (
                isA<calculatedFvPatchField<typename GeoField::value_type> >
                (pfld)
            )
            {
                tSubFld().boundaryField()[patchI] ==
                    pTraits<typename GeoField::value_type>::zero;
            }
        }

        // Store on subMesh
        GeoField* subFld = tSubFld.ptr();
        subFld->rename(fld.name());
        subFld->writeOpt() = IOobject::AUTO_WRITE;
        subFld->store();
    }
}


template<class GeoField>
void subsetSurfaceFields
(
    const fvMesh& mesh,
    const fvMesh& subMesh,
    const labelList& faceMap
)
{
    const labelList patchMap(identity(mesh.boundaryMesh().size()));

    HashTable<const GeoField*> fields
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );
    forAllConstIter(typename HashTable<const GeoField*>, fields, iter)
    {
        const GeoField& fld = *iter();

        Info<< "Mapping field " << fld.name() << endl;

        tmp<GeoField> tSubFld
        (
            fvMeshSubset::interpolate
            (
                fld,
                subMesh,
                patchMap,
                faceMap
            )
        );

        // Hack: set value to 0 for introduced patches (since don't
        //       get initialised.
        forAll(tSubFld().boundaryField(), patchI)
        {
            const fvsPatchField<typename GeoField::value_type>& pfld =
                tSubFld().boundaryField()[patchI];

            if
            (
                isA<calculatedFvsPatchField<typename GeoField::value_type> >
                (pfld)
            )
            {
                tSubFld().boundaryField()[patchI] ==
                    pTraits<typename GeoField::value_type>::zero;
            }
        }

        // Store on subMesh
        GeoField* subFld = tSubFld.ptr();
        subFld->rename(fld.name());
        subFld->writeOpt() = IOobject::AUTO_WRITE;
        subFld->store();
    }
}

// Select all cells not in the region
labelList getNonRegionCells(const labelList& cellRegion, const label regionI)
{
    DynamicList<label> nonRegionCells(cellRegion.size());
    forAll(cellRegion, cellI)
    {
        if (cellRegion[cellI] != regionI)
        {
            nonRegionCells.append(cellI);
        }
    }
    return nonRegionCells.shrink();
}


// Get per region-region interface the sizes. If sumParallel sums sizes.
// Returns interfaces as straight list for looping in identical order.
void getInterfaceSizes
(
    const polyMesh& mesh,
    const labelList& cellRegion,
    const bool sumParallel,

    edgeList& interfaces,
    EdgeMap<label>& interfaceSizes
)
{

    // Internal faces
    // ~~~~~~~~~~~~~~

    forAll(mesh.faceNeighbour(), faceI)
    {
        label ownRegion = cellRegion[mesh.faceOwner()[faceI]];
        label neiRegion = cellRegion[mesh.faceNeighbour()[faceI]];

        if (ownRegion != neiRegion)
        {
            edge interface
            (
                min(ownRegion, neiRegion),
                max(ownRegion, neiRegion)
            );

            EdgeMap<label>::iterator iter = interfaceSizes.find(interface);

            if (iter != interfaceSizes.end())
            {
                iter()++;
            }
            else
            {
                interfaceSizes.insert(interface, 1);
            }
        }
    }

    // Boundary faces
    // ~~~~~~~~~~~~~~

    // Neighbour cellRegion.
    labelList coupledRegion(mesh.nFaces()-mesh.nInternalFaces());

    forAll(coupledRegion, i)
    {
        label cellI = mesh.faceOwner()[i+mesh.nInternalFaces()];
        coupledRegion[i] = cellRegion[cellI];
    }
    syncTools::swapBoundaryFaceList(mesh, coupledRegion, false);

    forAll(coupledRegion, i)
    {
        label faceI = i+mesh.nInternalFaces();
        label ownRegion = cellRegion[mesh.faceOwner()[faceI]];
        label neiRegion = coupledRegion[i];

        if (ownRegion != neiRegion)
        {
            edge interface
            (
                min(ownRegion, neiRegion),
                max(ownRegion, neiRegion)
            );

            EdgeMap<label>::iterator iter = interfaceSizes.find(interface);

            if (iter != interfaceSizes.end())
            {
                iter()++;
            }
            else
            {
                interfaceSizes.insert(interface, 1);
            }
        }
    }


    if (sumParallel && Pstream::parRun())
    {
        if (Pstream::master())
        {
            // Receive and add to my sizes
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                IPstream fromSlave(Pstream::blocking, slave);

                EdgeMap<label> slaveSizes(fromSlave);

                forAllConstIter(EdgeMap<label>, slaveSizes, slaveIter)
                {
                    EdgeMap<label>::iterator masterIter =
                        interfaceSizes.find(slaveIter.key());

                    if (masterIter != interfaceSizes.end())
                    {
                        masterIter() += slaveIter();
                    }
                    else
                    {
                        interfaceSizes.insert(slaveIter.key(), slaveIter());
                    }
                }
            }

            // Distribute
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                // Receive the edges using shared points from the slave.
                OPstream toSlave(Pstream::blocking, slave);
                toSlave << interfaceSizes;
            }
        }
        else
        {
            // Send to master
            {
                OPstream toMaster(Pstream::blocking, Pstream::masterNo());
                toMaster << interfaceSizes;
            }
            // Receive from master
            {
                IPstream fromMaster(Pstream::blocking, Pstream::masterNo());
                fromMaster >> interfaceSizes;
            }
        }
    }

    // Make sure all processors have interfaces in same order
    interfaces = interfaceSizes.toc();
    if (sumParallel)
    {
        Pstream::scatter(interfaces);
    }
}


// Create mesh for region.
autoPtr<mapPolyMesh> createRegionMesh
(
    const labelList& cellRegion,
    const EdgeMap<label>& interfaceToPatch,
    const fvMesh& mesh,
    const label regionI,
    const word& regionName,
    autoPtr<fvMesh>& newMesh
)
{
    // Create dummy system/fv*
    {
        IOobject io
        (
            "fvSchemes",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        Info<< "Testing:" << io.objectPath() << endl;

        if (!io.headerOk())
        // if (!exists(io.objectPath()))
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            dictionary divDict;
            dummyDict.add("divSchemes", divDict);
            dictionary gradDict;
            dummyDict.add("gradSchemes", gradDict);
            dictionary laplDict;
            dummyDict.add("laplacianSchemes", laplDict);

            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }
    {
        IOobject io
        (
            "fvSolution",
            mesh.time().system(),
            regionName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!io.headerOk())
        //if (!exists(io.objectPath()))
        {
            Info<< "Writing dummy " << regionName/io.name() << endl;
            dictionary dummyDict;
            IOdictionary(io, dummyDict).regIOobject::write();
        }
    }


    // Neighbour cellRegion.
    labelList coupledRegion(mesh.nFaces()-mesh.nInternalFaces());

    forAll(coupledRegion, i)
    {
        label cellI = mesh.faceOwner()[i+mesh.nInternalFaces()];
        coupledRegion[i] = cellRegion[cellI];
    }
    syncTools::swapBoundaryFaceList(mesh, coupledRegion, false);


    // Topology change container. Start off from existing mesh.
    polyTopoChange meshMod(mesh);

    // Cell remover engine
    removeCells cellRemover(mesh);

    // Select all but region cells
    labelList cellsToRemove(getNonRegionCells(cellRegion, regionI));

    // Find out which faces will get exposed. Note that this
    // gets faces in mesh face order. So both regions will get same
    // face in same order (important!)
    labelList exposedFaces = cellRemover.getExposedFaces(cellsToRemove);

    labelList exposedPatchIDs(exposedFaces.size());
    forAll(exposedFaces, i)
    {
        label faceI = exposedFaces[i];

        label ownRegion = cellRegion[mesh.faceOwner()[faceI]];
        label neiRegion = -1;

        if (mesh.isInternalFace(faceI))
        {
            neiRegion = cellRegion[mesh.faceNeighbour()[faceI]];
        }
        else
        {
            neiRegion = coupledRegion[faceI-mesh.nInternalFaces()];
        }

        label otherRegion = -1;

        if (ownRegion == regionI && neiRegion != regionI)
        {
            otherRegion = neiRegion;
        }
        else if (ownRegion != regionI && neiRegion == regionI)
        {
            otherRegion = ownRegion;
        }
        else
        {
            FatalErrorIn("createRegionMesh(..)")
                << "Exposed face:" << faceI
                << " fc:" << mesh.faceCentres()[faceI]
                << " has owner region " << ownRegion
                << " and neighbour region " << neiRegion
                << " when handling region:" << regionI
                << exit(FatalError);
        }

        if (otherRegion != -1)
        {
            edge interface(regionI, otherRegion);

            // Find the patch.
            if (regionI < otherRegion)
            {
                exposedPatchIDs[i] = interfaceToPatch[interface];
            }
            else
            {
                exposedPatchIDs[i] = interfaceToPatch[interface]+1;
            }
        }
    }

    // Remove faces
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        meshMod
    );

    autoPtr<mapPolyMesh> map = meshMod.makeMesh
    (
        newMesh,
        IOobject
        (
            regionName,
            mesh.time().timeName(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    return map;
}


void createAndWriteRegion
(
    const fvMesh& mesh,
    const labelList& cellRegion,
    const wordList& regionNames,
    const EdgeMap<label>& interfaceToPatch,
    const label regionI,
    const word& newMeshInstance
)
{
    Info<< "Creating mesh for region " << regionI
        << ' ' << regionNames[regionI] << endl;

    autoPtr<fvMesh> newMesh;
    autoPtr<mapPolyMesh> map = createRegionMesh
    (
        cellRegion,
        interfaceToPatch,
        mesh,
        regionI,
        regionNames[regionI],
        newMesh
    );

    Info<< "Mapping fields" << endl;

    // Map existing fields
    newMesh().updateMesh(map());

    // Add subsetted fields
    subsetVolFields<volScalarField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap()
    );
    subsetVolFields<volVectorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap()
    );
    subsetVolFields<volSphericalTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap()
    );
    subsetVolFields<volSymmTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap()
    );
    subsetVolFields<volTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap()
    );

    subsetSurfaceFields<surfaceScalarField>
    (
        mesh,
        newMesh(),
        map().faceMap()
    );
    subsetSurfaceFields<surfaceVectorField>
    (
        mesh,
        newMesh(),
        map().faceMap()
    );
    subsetSurfaceFields<surfaceSphericalTensorField>
    (
        mesh,
        newMesh(),
        map().faceMap()
    );
    subsetSurfaceFields<surfaceSymmTensorField>
    (
        mesh,
        newMesh(),
        map().faceMap()
    );
    subsetSurfaceFields<surfaceTensorField>
    (
        mesh,
        newMesh(),
        map().faceMap()
    );


    const polyBoundaryMesh& newPatches = newMesh().boundaryMesh();
    newPatches.checkParallelSync(true);

    // Delete empty patches
    // ~~~~~~~~~~~~~~~~~~~~

    // Create reordering list to move patches-to-be-deleted to end
    labelList oldToNew(newPatches.size(), -1);
    label newI = 0;

    Info<< "Deleting empty patches" << endl;

    // Assumes all non-proc boundaries are on all processors!
    forAll(newPatches, patchI)
    {
        const polyPatch& pp = newPatches[patchI];

        if (!isA<processorPolyPatch>(pp))
        {
            if (returnReduce(pp.size(), sumOp<label>()) > 0)
            {
                oldToNew[patchI] = newI++;
            }
        }
    }

    // Same for processor patches (but need no reduction)
    forAll(newPatches, patchI)
    {
        const polyPatch& pp = newPatches[patchI];

        if (isA<processorPolyPatch>(pp) && pp.size())
        {
            oldToNew[patchI] = newI++;
        }
    }

    const label nNewPatches = newI;

    // Move all deleteable patches to the end
    forAll(oldToNew, patchI)
    {
        if (oldToNew[patchI] == -1)
        {
            oldToNew[patchI] = newI++;
        }
    }

    reorderPatches(newMesh(), oldToNew, nNewPatches);


    Info<< "Writing new mesh" << endl;

    newMesh().setInstance(newMeshInstance);
    newMesh().write();

    // Write addressing files like decomposePar
    Info<< "Writing addressing to base mesh" << endl;

    labelIOList pointProcAddressing
    (
        IOobject
        (
            "pointRegionAddressing",
            newMesh().facesInstance(),
            newMesh().meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        map().pointMap()
    );
    Info<< "Writing map " << pointProcAddressing.name()
        << " from region" << regionI
        << " points back to base mesh." << endl;
    pointProcAddressing.write();

    labelIOList faceProcAddressing
    (
        IOobject
        (
            "faceRegionAddressing",
            newMesh().facesInstance(),
            newMesh().meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        newMesh().nFaces()
    );
    forAll(faceProcAddressing, faceI)
    {
        // face + turning index. (see decomposePar)
        // Is the face pointing in the same direction?
        label oldFaceI = map().faceMap()[faceI];

        if
        (
            map().cellMap()[newMesh().faceOwner()[faceI]]
         == mesh.faceOwner()[oldFaceI]
        )
        {
            faceProcAddressing[faceI] = oldFaceI+1;
        }
        else
        {
            faceProcAddressing[faceI] = -(oldFaceI+1);
        }
    }
    Info<< "Writing map " << faceProcAddressing.name()
        << " from region" << regionI
        << " faces back to base mesh." << endl;
    faceProcAddressing.write();

    labelIOList cellProcAddressing
    (
        IOobject
        (
            "cellRegionAddressing",
            newMesh().facesInstance(),
            newMesh().meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        map().cellMap()
    );
    Info<< "Writing map " <<cellProcAddressing.name()
        << " from region" << regionI
        << " cells back to base mesh." << endl;
    cellProcAddressing.write();
}


// Create for every region-region interface with non-zero size two patches.
// First one is for minimumregion to maximumregion.
// Note that patches get created in same order on all processors (if parallel)
// since looping over synchronised 'interfaces'.
EdgeMap<label> addRegionPatches
(
    fvMesh& mesh,
    const labelList& cellRegion,
    const label nCellRegions,
    const edgeList& interfaces,
    const EdgeMap<label>& interfaceSizes,
    const wordList& regionNames
)
{
    // Check that all patches are present in same order.
    mesh.boundaryMesh().checkParallelSync(true);

    Info<< nl << "Adding patches" << nl << endl;

    EdgeMap<label> interfaceToPatch(nCellRegions);

    forAll(interfaces, interI)
    {
        const edge& e = interfaces[interI];

        if (interfaceSizes[e] > 0)
        {
            const word inter1 = regionNames[e[0]] + "_to_" + regionNames[e[1]];
            const word inter2 = regionNames[e[1]] + "_to_" + regionNames[e[0]];

            directMappedWallPolyPatch patch1
            (
                inter1,
                0,                  // overridden
                0,                  // overridden
                0,                  // overridden
                regionNames[e[1]],  // sampleRegion
                directMappedPatchBase::NEARESTPATCHFACE,
                inter2,             // samplePatch
                point::zero,        // offset
                mesh.boundaryMesh()
            );

            label patchI = addPatch(mesh, patch1);

            directMappedWallPolyPatch patch2
            (
                inter2,
                0,
                0,
                0,
                regionNames[e[0]],  // sampleRegion
                directMappedPatchBase::NEARESTPATCHFACE,
                inter1,
                point::zero,        // offset
                mesh.boundaryMesh()
            );
            addPatch(mesh, patch2);

            Info<< "For interface between region " << e[0]
                << " and " << e[1] << " added patch " << patchI
                << " " << mesh.boundaryMesh()[patchI].name()
                << endl;

            interfaceToPatch.insert(e, patchI);
        }
    }
    return interfaceToPatch;
}


// Find region that covers most of cell zone
label findCorrespondingRegion
(
    const labelList& existingZoneID,    // per cell the (unique) zoneID
    const labelList& cellRegion,
    const label nCellRegions,
    const label zoneI,
    const label minOverlapSize
)
{
    // Per region the number of cells in zoneI
    labelList cellsInZone(nCellRegions, 0);

    forAll(cellRegion, cellI)
    {
        if (existingZoneID[cellI] == zoneI)
        {
            cellsInZone[cellRegion[cellI]]++;
        }
    }

    Pstream::listCombineGather(cellsInZone, plusEqOp<label>());
    Pstream::listCombineScatter(cellsInZone);

    // Pick region with largest overlap of zoneI
    label regionI = findMax(cellsInZone);


    if (cellsInZone[regionI] < minOverlapSize)
    {
        // Region covers too little of zone. Not good enough.
        regionI = -1;
    }
    else
    {
        // Check that region contains no cells that aren't in cellZone.
        forAll(cellRegion, cellI)
        {
            if (cellRegion[cellI] == regionI && existingZoneID[cellI] != zoneI)
            {
                // cellI in regionI but not in zoneI
                regionI = -1;
                break;
            }
        }
        // If one in error, all should be in error. Note that branch gets taken
        // on all procs.
        reduce(regionI, minOp<label>());
    }

    return regionI;
}


//// Checks if cellZone has corresponding cellRegion.
//label findCorrespondingRegion
//(
//    const cellZoneMesh& cellZones,
//    const labelList& existingZoneID,    // per cell the (unique) zoneID
//    const labelList& cellRegion,
//    const label nCellRegions,
//    const label zoneI
//)
//{
//    // Region corresponding to zone. Start off with special value: no
//    // corresponding region.
//    label regionI = labelMax;
//
//    const cellZone& cz = cellZones[zoneI];
//
//    if (cz.empty())
//    {
//        // My local portion is empty. Maps to any empty cellZone. Mark with
//        // special value which can get overwritten by other processors.
//        regionI = -1;
//    }
//    else
//    {
//        regionI = cellRegion[cz[0]];
//
//        forAll(cz, i)
//        {
//            if (cellRegion[cz[i]] != regionI)
//            {
//                regionI = labelMax;
//                break;
//            }
//        }
//    }
//
//    // Determine same zone over all processors.
//    reduce(regionI, maxOp<label>());
//
//
//    // 2. All of region present?
//
//    if (regionI == labelMax)
//    {
//        regionI = -1;
//    }
//    else if (regionI != -1)
//    {
//        forAll(cellRegion, cellI)
//        {
//            if (cellRegion[cellI] == regionI && existingZoneID[cellI] != zoneI)
//            {
//                // cellI in regionI but not in zoneI
//                regionI = -1;
//                break;
//            }
//        }
//        // If one in error, all should be in error. Note that branch gets taken
//        // on all procs.
//        reduce(regionI, minOp<label>());
//    }
//
//    return regionI;
//}


// Get zone per cell
// - non-unique zoning
// - coupled zones
void getZoneID
(
    const polyMesh& mesh,
    const cellZoneMesh& cellZones,
    labelList& zoneID,
    labelList& neiZoneID
)
{
    // Existing zoneID
    zoneID.setSize(mesh.nCells());
    zoneID = -1;

    forAll(cellZones, zoneI)
    {
        const cellZone& cz = cellZones[zoneI];

        forAll(cz, i)
        {
            label cellI = cz[i];
            if (zoneID[cellI] == -1)
            {
                zoneID[cellI] = zoneI;
            }
            else
            {
                FatalErrorIn("getZoneID(..)")
                    << "Cell " << cellI << " with cell centre "
                    << mesh.cellCentres()[cellI]
                    << " is multiple zones. This is not allowed." << endl
                    << "It is in zone " << cellZones[zoneID[cellI]].name()
                    << " and in zone " << cellZones[zoneI].name()
                    << exit(FatalError);
            }
        }
    }

    // Neighbour zoneID.
    neiZoneID.setSize(mesh.nFaces()-mesh.nInternalFaces());

    forAll(neiZoneID, i)
    {
        neiZoneID[i] = zoneID[mesh.faceOwner()[i+mesh.nInternalFaces()]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiZoneID, false);
}


void matchRegions
(
    const bool sloppyCellZones,
    const polyMesh& mesh,

    const label nCellRegions,
    const labelList& cellRegion,

    labelList& regionToZone,
    wordList& regionNames,
    labelList& zoneToRegion
)
{
    const cellZoneMesh& cellZones = mesh.cellZones();

    regionToZone.setSize(nCellRegions, -1);
    regionNames.setSize(nCellRegions);
    zoneToRegion.setSize(cellZones.size(), -1);

    // Get current per cell zoneID
    labelList zoneID(mesh.nCells(), -1);
    labelList neiZoneID(mesh.nFaces()-mesh.nInternalFaces());
    getZoneID(mesh, cellZones, zoneID, neiZoneID);

    // Sizes per cellzone
    labelList zoneSizes(cellZones.size(), 0);
    {
        List<wordList> zoneNames(Pstream::nProcs());
        zoneNames[Pstream::myProcNo()] = cellZones.names();
        Pstream::gatherList(zoneNames);
        Pstream::scatterList(zoneNames);

        forAll(zoneNames, procI)
        {
            if (zoneNames[procI] != zoneNames[0])
            {
                FatalErrorIn("matchRegions(..)")
                    << "cellZones not synchronised across processors." << endl
                    << "Master has cellZones " << zoneNames[0] << endl
                    << "Processor " << procI
                    << " has cellZones " << zoneNames[procI]
                    << exit(FatalError);
            }
        }

        forAll(cellZones, zoneI)
        {
            zoneSizes[zoneI] = returnReduce
            (
                cellZones[zoneI].size(),
                sumOp<label>()
            );
        }
    }


    if (sloppyCellZones)
    {
        Info<< "Trying to match regions to existing cell zones;"
            << " region can be subset of cell zone." << nl << endl;

        forAll(cellZones, zoneI)
        {
            label regionI = findCorrespondingRegion
            (
                zoneID,
                cellRegion,
                nCellRegions,
                zoneI,
                label(0.5*zoneSizes[zoneI]) // minimum overlap
            );

            if (regionI != -1)
            {
                Info<< "Sloppily matched region " << regionI
                    //<< " size " << regionSizes[regionI]
                    << " to zone " << zoneI << " size " << zoneSizes[zoneI]
                    << endl;
                zoneToRegion[zoneI] = regionI;
                regionToZone[regionI] = zoneI;
                regionNames[regionI] = cellZones[zoneI].name();
            }
        }
    }
    else
    {
        Info<< "Trying to match regions to existing cell zones." << nl << endl;

        forAll(cellZones, zoneI)
        {
            label regionI = findCorrespondingRegion
            (
                zoneID,
                cellRegion,
                nCellRegions,
                zoneI,
                1               // minimum overlap
            );

            if (regionI != -1)
            {
                zoneToRegion[zoneI] = regionI;
                regionToZone[regionI] = zoneI;
                regionNames[regionI] = cellZones[zoneI].name();
            }
        }
    }
    // Allocate region names for unmatched regions.
    forAll(regionToZone, regionI)
    {
        if (regionToZone[regionI] == -1)
        {
            regionNames[regionI] = "domain" + Foam::name(regionI);
        }
    }
}


void writeCellToRegion(const fvMesh& mesh, const labelList& cellRegion)
{
    // Write to manual decomposition option
    {
        labelIOList cellToRegion
        (
            IOobject
            (
                "cellToRegion",
                mesh.facesInstance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            cellRegion
        );
        cellToRegion.write();

        Info<< "Writing region per cell file (for manual decomposition) to "
            << cellToRegion.objectPath() << nl << endl;
    }
    // Write for postprocessing
    {
        volScalarField cellToRegion
        (
            IOobject
            (
                "cellToRegion",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );
        forAll(cellRegion, cellI)
        {
            cellToRegion[cellI] = cellRegion[cellI];
        }
        cellToRegion.write();

        Info<< "Writing region per cell as volScalarField to "
            << cellToRegion.objectPath() << nl << endl;
    }
}



// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("cellZones", "");
    argList::validOptions.insert("cellZonesOnly", "");
    argList::validOptions.insert("cellZonesFileOnly", "cellZonesName");
    argList::validOptions.insert("blockedFaces", "faceSet");
    argList::validOptions.insert("makeCellZones", "");
    argList::validOptions.insert("largestOnly", "");
    argList::validOptions.insert("insidePoint", "point");
    argList::validOptions.insert("overwrite", "");
    argList::validOptions.insert("detectOnly", "");
    argList::validOptions.insert("sloppyCellZones", "");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createMesh.H"
    const word oldInstance = mesh.pointsInstance();

    word blockedFacesName;
    if (args.optionFound("blockedFaces"))
    {
        blockedFacesName = args.option("blockedFaces");
        Info<< "Reading blocked internal faces from faceSet "
            << blockedFacesName << nl << endl;
    }

    bool makeCellZones    = args.optionFound("makeCellZones");
    bool largestOnly      = args.optionFound("largestOnly");
    bool insidePoint      = args.optionFound("insidePoint");
    bool useCellZones     = args.optionFound("cellZones");
    bool useCellZonesOnly = args.optionFound("cellZonesOnly");
    bool useCellZonesFile = args.optionFound("cellZonesFileOnly");
    bool overwrite        = args.optionFound("overwrite");
    bool detectOnly       = args.optionFound("detectOnly");
    bool sloppyCellZones  = args.optionFound("sloppyCellZones");

    if
    (
        (useCellZonesOnly || useCellZonesFile)
     && (
            blockedFacesName != word::null
         || useCellZones
        )
    )
    {
        FatalErrorIn(args.executable())
            << "You cannot specify both -cellZonesOnly or -cellZonesFileOnly"
            << " (which specify complete split)"
            << " in combination with -blockedFaces or -cellZones"
            << " (which imply a split based on topology)"
            << exit(FatalError);
    }



    if (insidePoint && largestOnly)
    {
        FatalErrorIn(args.executable())
            << "You cannot specify both -largestOnly"
            << " (keep region with most cells)"
            << " and -insidePoint (keep region containing point)"
            << exit(FatalError);
    }


    const cellZoneMesh& cellZones = mesh.cellZones();

    // Existing zoneID
    labelList zoneID(mesh.nCells(), -1);
    // Neighbour zoneID.
    labelList neiZoneID(mesh.nFaces()-mesh.nInternalFaces());
    getZoneID(mesh, cellZones, zoneID, neiZoneID);



    // Determine per cell the region it belongs to
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // cellRegion is the labelList with the region per cell.
    labelList cellRegion;
    // Region per zone
    labelList regionToZone;
    // Name of region
    wordList regionNames;
    // Zone to region
    labelList zoneToRegion;

    label nCellRegions = 0;
    if (useCellZonesOnly)
    {
        Info<< "Using current cellZones to split mesh into regions."
            << " This requires all"
            << " cells to be in one and only one cellZone." << nl << endl;

        label unzonedCellI = findIndex(zoneID, -1);
        if (unzonedCellI != -1)
        {
            FatalErrorIn(args.executable())
                << "For the cellZonesOnly option all cells "
                << "have to be in a cellZone." << endl
                << "Cell " << unzonedCellI
                << " at" << mesh.cellCentres()[unzonedCellI]
                << " is not in a cellZone. There might be more unzoned cells."
                << exit(FatalError);
        }
        cellRegion = zoneID;
        nCellRegions = gMax(cellRegion)+1;
        regionToZone.setSize(nCellRegions);
        regionNames.setSize(nCellRegions);
        zoneToRegion.setSize(cellZones.size(), -1);
        for (label regionI = 0; regionI < nCellRegions; regionI++)
        {
            regionToZone[regionI] = regionI;
            zoneToRegion[regionI] = regionI;
            regionNames[regionI] = cellZones[regionI].name();
        }
    }
    else if (useCellZonesFile)
    {
        const word zoneFile = args.option("cellZonesFileOnly");
        Info<< "Reading split from cellZones file " << zoneFile << endl
            << "This requires all"
            << " cells to be in one and only one cellZone." << nl << endl;

        cellZoneMesh newCellZones
        (
            IOobject
            (
                zoneFile,
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh
        );

        labelList newZoneID(mesh.nCells(), -1);
        labelList newNeiZoneID(mesh.nFaces()-mesh.nInternalFaces());
        getZoneID(mesh, newCellZones, newZoneID, newNeiZoneID);

        label unzonedCellI = findIndex(newZoneID, -1);
        if (unzonedCellI != -1)
        {
            FatalErrorIn(args.executable())
                << "For the cellZonesFileOnly option all cells "
                << "have to be in a cellZone." << endl
                << "Cell " << unzonedCellI
                << " at" << mesh.cellCentres()[unzonedCellI]
                << " is not in a cellZone. There might be more unzoned cells."
                << exit(FatalError);
        }
        cellRegion = newZoneID;
        nCellRegions = gMax(cellRegion)+1;
        zoneToRegion.setSize(newCellZones.size(), -1);
        regionToZone.setSize(nCellRegions);
        regionNames.setSize(nCellRegions);
        for (label regionI = 0; regionI < nCellRegions; regionI++)
        {
            regionToZone[regionI] = regionI;
            zoneToRegion[regionI] = regionI;
            regionNames[regionI] = newCellZones[regionI].name();
        }
    }
    else
    {
        // Determine connected regions
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Mark additional faces that are blocked
        boolList blockedFace;

        // Read from faceSet
        if (blockedFacesName.size())
        {
            faceSet blockedFaceSet(mesh, blockedFacesName);
            Info<< "Read "
                << returnReduce(blockedFaceSet.size(), sumOp<label>())
                << " blocked faces from set " << blockedFacesName << nl << endl;

            blockedFace.setSize(mesh.nFaces(), false);

            forAllConstIter(faceSet, blockedFaceSet, iter)
            {
                blockedFace[iter.key()] = true;
            }
        }

        // Imply from differing cellZones
        if (useCellZones)
        {
            blockedFace.setSize(mesh.nFaces(), false);

            for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
            {
                label own = mesh.faceOwner()[faceI];
                label nei = mesh.faceNeighbour()[faceI];

                if (zoneID[own] != zoneID[nei])
                {
                    blockedFace[faceI] = true;
                }
            }

            // Different cellZones on either side of processor patch.
            forAll(neiZoneID, i)
            {
                label faceI = i+mesh.nInternalFaces();

                if (zoneID[mesh.faceOwner()[faceI]] != neiZoneID[i])
                {
                    blockedFace[faceI] = true;
                }
            }
        }

        // Do a topological walk to determine regions
        regionSplit regions(mesh, blockedFace);
        nCellRegions = regions.nRegions();
        cellRegion.transfer(regions);

        // Make up region names. If possible match them to existing zones.
        matchRegions
        (
            sloppyCellZones,
            mesh,
            nCellRegions,
            cellRegion,

            regionToZone,
            regionNames,
            zoneToRegion
        );
    }

    Info<< endl << "Number of regions:" << nCellRegions << nl << endl;


    // Write decomposition to file
    writeCellToRegion(mesh, cellRegion);



    // Sizes per region
    // ~~~~~~~~~~~~~~~~

    labelList regionSizes(nCellRegions, 0);

    forAll(cellRegion, cellI)
    {
        regionSizes[cellRegion[cellI]]++;
    }
    forAll(regionSizes, regionI)
    {
        reduce(regionSizes[regionI], sumOp<label>());
    }

    Info<< "Region\tCells" << nl
        << "------\t-----" << endl;

    forAll(regionSizes, regionI)
    {
        Info<< regionI << '\t' << regionSizes[regionI] << nl;
    }
    Info<< endl;



    // Print region to zone
    Info<< "Region\tZone\tName" << nl
        << "------\t----\t----" << endl;
    forAll(regionToZone, regionI)
    {
        Info<< regionI << '\t' << regionToZone[regionI] << '\t'
            << regionNames[regionI] << nl;
    }
    Info<< endl;



    // Since we're going to mess with patches make sure all non-processor ones
    // are on all processors.
    mesh.boundaryMesh().checkParallelSync(true);


    // Sizes of interface between regions. From pair of regions to number of
    // faces.
    edgeList interfaces;
    EdgeMap<label> interfaceSizes;
    getInterfaceSizes
    (
        mesh,
        cellRegion,
        true,      // sum in parallel?

        interfaces,
        interfaceSizes
    );

    Info<< "Sizes inbetween regions:" << nl << nl
        << "Region\tRegion\tFaces" << nl
        << "------\t------\t-----" << endl;

    forAll(interfaces, interI)
    {
        const edge& e = interfaces[interI];

        Info<< e[0] << '\t' << e[1] << '\t' << interfaceSizes[e] << nl;
    }
    Info<< endl;


    if (detectOnly)
    {
        return 0;
    }


    // Read objects in time directory
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    ReadFields(mesh, objects, stFlds);

    Info<< endl;


    // Remove any demand-driven fields ('S', 'V' etc)
    mesh.clearOut();


    if (nCellRegions == 1)
    {
        Info<< "Only one region. Doing nothing." << endl;
    }
    else if (makeCellZones)
    {
        Info<< "Putting cells into cellZones instead of splitting mesh."
            << endl;

        // Check if region overlaps with existing zone. If so keep.

        for (label regionI = 0; regionI < nCellRegions; regionI++)
        {
            label zoneI = regionToZone[regionI];

            if (zoneI != -1)
            {
                Info<< "    Region " << regionI << " : corresponds to existing"
                    << " cellZone "
                    << zoneI << ' ' << cellZones[zoneI].name() << endl;
            }
            else
            {
                // Create new cellZone.
                labelList regionCells = findIndices(cellRegion, regionI);

                word zoneName = "region" + Foam::name(regionI);

                zoneI = cellZones.findZoneID(zoneName);

                if (zoneI == -1)
                {
                    zoneI = cellZones.size();
                    mesh.cellZones().setSize(zoneI+1);
                    mesh.cellZones().set
                    (
                        zoneI,
                        new cellZone
                        (
                            zoneName,           //name
                            regionCells,        //addressing
                            zoneI,              //index
                            cellZones           //cellZoneMesh
                        )
                    );
                }
                else
                {
                    mesh.cellZones()[zoneI].clearAddressing();
                    mesh.cellZones()[zoneI] = regionCells;
                }
                Info<< "    Region " << regionI << " : created new cellZone "
                    << zoneI << ' ' << cellZones[zoneI].name() << endl;
            }
        }
        mesh.cellZones().writeOpt() = IOobject::AUTO_WRITE;

        if (!overwrite)
        {
            runTime++;
            mesh.setInstance(runTime.timeName());
        }
        else
        {
            mesh.setInstance(oldInstance);
        }

        Info<< "Writing cellZones as new mesh to time " << runTime.timeName()
            << nl << endl;

        mesh.write();


        // Write cellSets for convenience
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Info<< "Writing cellSets corresponding to cellZones." << nl << endl;

        forAll(cellZones, zoneI)
        {
            const cellZone& cz = cellZones[zoneI];

            cellSet(mesh, cz.name(), cz).write();
        }
    }
    else
    {
        // Add patches for interfaces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Add all possible patches. Empty ones get filtered later on.
        Info<< nl << "Adding patches" << nl << endl;

        EdgeMap<label> interfaceToPatch
        (
            addRegionPatches
            (
                mesh,
                cellRegion,
                nCellRegions,
                interfaces,
                interfaceSizes,
                regionNames
            )
        );


        if (!overwrite)
        {
            runTime++;
        }


        // Create regions
        // ~~~~~~~~~~~~~~

        if (insidePoint)
        {
            point insidePoint(args.optionLookup("insidePoint")());

            label regionI = -1;

            label cellI = mesh.findCell(insidePoint);

            Info<< nl << "Found point " << insidePoint << " in cell " << cellI
                << endl;

            if (cellI != -1)
            {
                regionI = cellRegion[cellI];
            }

            reduce(regionI, maxOp<label>());

            Info<< nl
                << "Subsetting region " << regionI
                << " containing point " << insidePoint << endl;

            if (regionI == -1)
            {
                FatalErrorIn(args.executable())
                    << "Point " << insidePoint
                    << " is not inside the mesh." << nl
                    << "Bounding box of the mesh:" << mesh.bounds()
                    << exit(FatalError);
            }

            createAndWriteRegion
            (
                mesh,
                cellRegion,
                regionNames,
                interfaceToPatch,
                regionI,
                (overwrite ? oldInstance : runTime.timeName())
            );
        }
        else if (largestOnly)
        {
            label regionI = findMax(regionSizes);

            Info<< nl
                << "Subsetting region " << regionI
                << " of size " << regionSizes[regionI] << endl;

            createAndWriteRegion
            (
                mesh,
                cellRegion,
                regionNames,
                interfaceToPatch,
                regionI,
                (overwrite ? oldInstance : runTime.timeName())
            );
        }
        else
        {
            // Split all
            for (label regionI = 0; regionI < nCellRegions; regionI++)
            {
                Info<< nl
                    << "Region " << regionI << nl
                    << "-------- " << endl;

                createAndWriteRegion
                (
                    mesh,
                    cellRegion,
                    regionNames,
                    interfaceToPatch,
                    regionI,
                    (overwrite ? oldInstance : runTime.timeName())
                );
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
