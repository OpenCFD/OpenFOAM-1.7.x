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

\*----------------------------------------------------------------------------*/

#include "CompactListList_dev.H"
#include "fvMeshDistribute.H"
#include "PstreamCombineReduceOps.H"
#include "fvMeshAdder.H"
#include "faceCoupleInfo.H"
#include "processorFvPatchField.H"
#include "processorFvsPatchField.H"
#include "polyTopoChange.H"
#include "removeCells.H"
#include "polyModifyFace.H"
#include "polyRemovePoint.H"
#include "mergePoints.H"
#include "mapDistributePolyMesh.H"
#include "surfaceFields.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::fvMeshDistribute, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::fvMeshDistribute::select
(
    const bool selectEqual,
    const labelList& values,
    const label value
)
{
    label n = 0;

    forAll(values, i)
    {
        if (selectEqual == (values[i] == value))
        {
            n++;
        }
    }

    labelList indices(n);
    n = 0;

    forAll(values, i)
    {
        if (selectEqual == (values[i] == value))
        {
            indices[n++] = i;
        }
    }
    return indices;
}


// Check all procs have same names and in exactly same order.
void Foam::fvMeshDistribute::checkEqualWordList
(
    const string& msg,
    const wordList& lst
)
{
    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = lst;
    Pstream::gatherList(allNames);
    Pstream::scatterList(allNames);

    for (label procI = 1; procI < Pstream::nProcs(); procI++)
    {
        if (allNames[procI] != allNames[0])
        {
            FatalErrorIn("fvMeshDistribute::checkEqualWordList(..)")
                << "When checking for equal " << msg.c_str() << " :" << endl
                << "processor0 has:" << allNames[0] << endl
                << "processor" << procI << " has:" << allNames[procI] << endl
                << msg.c_str() << " need to be synchronised on all processors."
                << exit(FatalError);
        }
    }
}


Foam::wordList Foam::fvMeshDistribute::mergeWordList(const wordList& procNames)
{
    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = procNames;
    Pstream::gatherList(allNames);
    Pstream::scatterList(allNames);

    HashSet<word> mergedNames;
    forAll(allNames, procI)
    {
        forAll(allNames[procI], i)
        {
            mergedNames.insert(allNames[procI][i]);
        }
    }
    return mergedNames.toc();
}


// Print some info on mesh.
void Foam::fvMeshDistribute::printMeshInfo(const fvMesh& mesh)
{
    Pout<< "Primitives:" << nl
        << "    points       :" << mesh.nPoints() << nl
        << "    internalFaces:" << mesh.nInternalFaces() << nl
        << "    faces        :" << mesh.nFaces() << nl
        << "    cells        :" << mesh.nCells() << nl;

    const fvBoundaryMesh& patches = mesh.boundary();

    Pout<< "Patches:" << endl;
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI].patch();

        Pout<< "    " << patchI << " name:" << pp.name()
            << " size:" << pp.size()
            << " start:" << pp.start()
            << " type:" << pp.type()
            << endl;
    }

    if (mesh.pointZones().size())
    {
        Pout<< "PointZones:" << endl;
        forAll(mesh.pointZones(), zoneI)
        {
            const pointZone& pz = mesh.pointZones()[zoneI];
            Pout<< "    " << zoneI << " name:" << pz.name()
                << " size:" << pz.size()
                << endl;
        }
    }
    if (mesh.faceZones().size())
    {
        Pout<< "FaceZones:" << endl;
        forAll(mesh.faceZones(), zoneI)
        {
            const faceZone& fz = mesh.faceZones()[zoneI];
            Pout<< "    " << zoneI << " name:" << fz.name()
                << " size:" << fz.size()
                << endl;
        }
    }
    if (mesh.cellZones().size())
    {
        Pout<< "CellZones:" << endl;
        forAll(mesh.cellZones(), zoneI)
        {
            const cellZone& cz = mesh.cellZones()[zoneI];
            Pout<< "    " << zoneI << " name:" << cz.name()
                << " size:" << cz.size()
                << endl;
        }
    }
}


void Foam::fvMeshDistribute::printCoupleInfo
(
    const primitiveMesh& mesh,
    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourceNewProc
)
{
    Pout<< nl
        << "Current coupling info:"
        << endl;

    forAll(sourceFace, bFaceI)
    {
        label meshFaceI = mesh.nInternalFaces() + bFaceI;

        Pout<< "    meshFace:" << meshFaceI
            << " fc:" << mesh.faceCentres()[meshFaceI]
            << " connects to proc:" << sourceProc[bFaceI]
            << "/face:" << sourceFace[bFaceI]
            << " which will move to proc:" << sourceNewProc[bFaceI]
            << endl;
    }
}


// Finds (non-empty) patch that exposed internal and proc faces can be put into.
Foam::label Foam::fvMeshDistribute::findNonEmptyPatch() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nonEmptyPatchI = -1;

    forAllReverse(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (!isA<emptyPolyPatch>(pp) && !pp.coupled())
        {
            nonEmptyPatchI = patchI;
            break;
        }
    }

    if (nonEmptyPatchI == -1)
    {
        FatalErrorIn("fvMeshDistribute::findNonEmptyPatch() const")
            << "Cannot find a patch which is neither of type empty nor"
            << " coupled in patches " << patches.names() << endl
            << "There has to be at least one such patch for"
            << " distribution to work" << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "findNonEmptyPatch : using patch " << nonEmptyPatchI
            << " name:" << patches[nonEmptyPatchI].name()
            << " type:" << patches[nonEmptyPatchI].type()
            << " to put exposed faces into." << endl;
    }


    // Do additional test for processor patches intermingled with non-proc
    // patches.
    label procPatchI = -1;

    forAll(patches, patchI)
    {
        if (isA<processorPolyPatch>(patches[patchI]))
        {
            procPatchI = patchI;
        }
        else if (procPatchI != -1)
        {
            FatalErrorIn("fvMeshDistribute::findNonEmptyPatch() const")
                << "Processor patches should be at end of patch list."
                << endl
                << "Have processor patch " << procPatchI
                << " followed by non-processor patch " << patchI
                << " in patches " << patches.names()
                << abort(FatalError);
        }
    }

    return nonEmptyPatchI;
}


// Appends processorPolyPatch. Returns patchID.
Foam::label Foam::fvMeshDistribute::addProcPatch
(
    const word& patchName,
    const label nbrProc
)
{
    // Clear local fields and e.g. polyMesh globalMeshData.
    mesh_.clearOut();


    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh_.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh_.boundary());

    if (polyPatches.findPatchID(patchName) != -1)
    {
        FatalErrorIn("fvMeshDistribute::addProcPatch(const word&, const label)")
            << "Cannot create patch " << patchName << " since already exists."
            << nl
            << "Current patch names:" << polyPatches.names()
            << exit(FatalError);
    }



    // Add the patch
    // ~~~~~~~~~~~~~

    label sz = polyPatches.size();

    // Add polyPatch
    polyPatches.setSize(sz+1);
    polyPatches.set
    (
        sz,
        new processorPolyPatch
        (
            patchName,
            0,              // size
            mesh_.nFaces(),
            sz,
            mesh_.boundaryMesh(),
            Pstream::myProcNo(),
            nbrProc
        )
    );
    fvPatches.setSize(sz+1);
    fvPatches.set
    (
        sz,
        fvPatch::New
        (
            polyPatches[sz],  // point to newly added polyPatch
            mesh_.boundary()
        )
    );

    return sz;
}


// Deletes last patch
void Foam::fvMeshDistribute::deleteTrailingPatch()
{
    // Clear local fields and e.g. polyMesh globalMeshData.
    mesh_.clearOut();

    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh_.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh_.boundary());

    if (polyPatches.empty())
    {
        FatalErrorIn("fvMeshDistribute::deleteTrailingPatch(fvMesh&)")
            << "No patches in mesh"
            << abort(FatalError);
    }

    label sz = polyPatches.size();

    label nFaces = polyPatches[sz-1].size();

    // Remove last polyPatch
    if (debug)
    {
        Pout<< "deleteTrailingPatch : Removing patch " << sz-1
            << " : " << polyPatches[sz-1].name() << " size:" << nFaces << endl;
    }

    if (nFaces)
    {
        FatalErrorIn("fvMeshDistribute::deleteTrailingPatch()")
            << "There are still " << nFaces << " faces in patch to be deleted "
            << sz-1 << ' ' << polyPatches[sz-1].name()
            << abort(FatalError);
    }

    // Remove actual patch
    polyPatches.setSize(sz-1);
    fvPatches.setSize(sz-1);
}


// Delete all processor patches. Move any processor faces into the last
// non-processor patch.
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::deleteProcPatches
(
    const label destinationPatch
)
{
    // New patchID per boundary faces to be repatched. Is -1 (no change)
    // or new patchID
    labelList newPatchID(mesh_.nFaces() - mesh_.nInternalFaces(), -1);

    label nProcPatches = 0;

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            if (debug)
            {
                Pout<< "Moving all faces of patch " << pp.name()
                    << " into patch " << destinationPatch
                    << endl;
            }

            label offset = pp.start() - mesh_.nInternalFaces();

            forAll(pp, i)
            {
                newPatchID[offset+i] = destinationPatch;
            }

            nProcPatches++;
        }
    }

    // Note: order of boundary faces been kept the same since the
    // destinationPatch is at the end and we have visited the patches in
    // incremental order.
    labelListList dummyFaceMaps;
    autoPtr<mapPolyMesh> map = repatch(newPatchID, dummyFaceMaps);


    // Delete (now empty) processor patches.
    forAllReverse(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            deleteTrailingPatch();
            deleteTrailingPatchFields<volScalarField>();
            deleteTrailingPatchFields<volVectorField>();
            deleteTrailingPatchFields<volSphericalTensorField>();
            deleteTrailingPatchFields<volSymmTensorField>();
            deleteTrailingPatchFields<volTensorField>();

            deleteTrailingPatchFields<surfaceScalarField>();
            deleteTrailingPatchFields<surfaceVectorField>();
            deleteTrailingPatchFields<surfaceSphericalTensorField>();
            deleteTrailingPatchFields<surfaceSymmTensorField>();
            deleteTrailingPatchFields<surfaceTensorField>();
        }
    }

    return map;
}


// Repatch the mesh.
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::repatch
(
    const labelList& newPatchID,         // per boundary face -1 or new patchID
    labelListList& constructFaceMap
)
{
    polyTopoChange meshMod(mesh_);

    forAll(newPatchID, bFaceI)
    {
        if (newPatchID[bFaceI] != -1)
        {
            label faceI = mesh_.nInternalFaces() + bFaceI;

            label zoneID = mesh_.faceZones().whichZone(faceI);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = mesh_.faceZones()[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh_.faces()[faceI],       // modified face
                    faceI,                      // label of face
                    mesh_.faceOwner()[faceI],   // owner
                    -1,                         // neighbour
                    false,                      // face flip
                    newPatchID[bFaceI],         // patch for face
                    false,                      // remove from zone
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                )
            );
        }
    }


    // Do mapping of fields from one patchField to the other ourselves since
    // is currently not supported by updateMesh.

    // Store boundary fields (we only do this for surfaceFields)
    PtrList<FieldField<fvsPatchField, scalar> > sFlds;
    saveBoundaryFields<scalar, surfaceMesh>(sFlds);
    PtrList<FieldField<fvsPatchField, vector> > vFlds;
    saveBoundaryFields<vector, surfaceMesh>(vFlds);
    PtrList<FieldField<fvsPatchField, sphericalTensor> > sptFlds;
    saveBoundaryFields<sphericalTensor, surfaceMesh>(sptFlds);
    PtrList<FieldField<fvsPatchField, symmTensor> > sytFlds;
    saveBoundaryFields<symmTensor, surfaceMesh>(sytFlds);
    PtrList<FieldField<fvsPatchField, tensor> > tFlds;
    saveBoundaryFields<tensor, surfaceMesh>(tFlds);

    // Change the mesh (no inflation). Note: parallel comms allowed.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields. No inflation, parallel sync.
    mesh_.updateMesh(map);

    // Map patch fields using stored boundary fields. Note: assumes order
    // of fields has not changed in object registry!
    mapBoundaryFields<scalar, surfaceMesh>(map, sFlds);
    mapBoundaryFields<vector, surfaceMesh>(map, vFlds);
    mapBoundaryFields<sphericalTensor, surfaceMesh>(map, sptFlds);
    mapBoundaryFields<symmTensor, surfaceMesh>(map, sytFlds);
    mapBoundaryFields<tensor, surfaceMesh>(map, tFlds);


    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }

    // Adapt constructMaps.

    if (debug)
    {
        label index = findIndex(map().reverseFaceMap(), -1);

        if (index != -1)
        {
            FatalErrorIn
            (
                "fvMeshDistribute::repatch(const labelList&, labelListList&)"
            )   << "reverseFaceMap contains -1 at index:"
                << index << endl
                << "This means that the repatch operation was not just"
                << " a shuffle?" << abort(FatalError);
        }
    }

    forAll(constructFaceMap, procI)
    {
        inplaceRenumber(map().reverseFaceMap(), constructFaceMap[procI]);
    }


    return map;
}


// Detect shared points. Need processor patches to be present.
// Background: when adding bits of mesh one can get points which
// share the same position but are only detectable to be topologically
// the same point when doing parallel analysis. This routine will
// merge those points.
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::mergeSharedPoints
(
    labelListList& constructPointMap
)
{
    // Find out which sets of points get merged and create a map from
    // mesh point to unique point.
    Map<label> pointToMaster
    (
        fvMeshAdder::findSharedPoints
        (
            mesh_,
            mergeTol_
        )
    );

    if (returnReduce(pointToMaster.size(), sumOp<label>()) == 0)
    {
        return autoPtr<mapPolyMesh>(NULL);
    }

    polyTopoChange meshMod(mesh_);

    fvMeshAdder::mergePoints(mesh_, pointToMaster, meshMod);

    // Change the mesh (no inflation). Note: parallel comms allowed.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields. No inflation, parallel sync.
    mesh_.updateMesh(map);

    // Adapt constructMaps for merged points.
    forAll(constructPointMap, procI)
    {
        labelList& constructMap = constructPointMap[procI];

        forAll(constructMap, i)
        {
            label oldPointI = constructMap[i];

            label newPointI = map().reversePointMap()[oldPointI];

            if (newPointI < -1)
            {
                constructMap[i] = -newPointI-2;
            }
            else if (newPointI >= 0)
            {
                constructMap[i] = newPointI;
            }
            else
            {
                FatalErrorIn("fvMeshDistribute::mergeSharedPoints()")
                    << "Problem. oldPointI:" << oldPointI
                    << " newPointI:" << newPointI << abort(FatalError);
            }
        }
    }
    return map;
}


// Construct the local environment of all boundary faces.
void Foam::fvMeshDistribute::getNeighbourData
(
    const labelList& distribution,
    labelList& sourceFace,
    labelList& sourceProc,
    labelList& sourceNewProc
) const
{
    label nBnd = mesh_.nFaces() - mesh_.nInternalFaces();
    sourceFace.setSize(nBnd);
    sourceProc.setSize(nBnd);
    sourceNewProc.setSize(nBnd);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Get neighbouring meshFace labels and new processor of coupled boundaries.
    labelList nbrFaces(nBnd, -1);
    labelList nbrNewProc(nBnd, -1);

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            label offset = pp.start() - mesh_.nInternalFaces();

            // Mesh labels of faces on this side
            forAll(pp, i)
            {
                label bndI = offset + i;
                nbrFaces[bndI] = pp.start()+i;
            }

            // Which processor they will end up on
            SubList<label>(nbrNewProc, pp.size(), offset).assign
            (
                UIndirectList<label>(distribution, pp.faceCells())()
            );
        }
    }


    // Exchange the boundary data
    syncTools::swapBoundaryFaceList(mesh_, nbrFaces, false);
    syncTools::swapBoundaryFaceList(mesh_, nbrNewProc, false);


    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        label offset = pp.start() - mesh_.nInternalFaces();

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            // Check which of the two faces we store.

            if (Pstream::myProcNo() < procPatch.neighbProcNo())
            {
                // Use my local face labels
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = pp.start()+i;
                    sourceProc[bndI] = Pstream::myProcNo();
                    sourceNewProc[bndI] = nbrNewProc[bndI];
                }
            }
            else
            {
                // Use my neighbours face labels
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = nbrFaces[bndI];
                    sourceProc[bndI] = procPatch.neighbProcNo();
                    sourceNewProc[bndI] = nbrNewProc[bndI];
                }
            }
        }
        else
        {
            // Normal (physical) boundary
            forAll(pp, i)
            {
                label bndI = offset + i;
                sourceFace[bndI] = patchI;
                sourceProc[bndI] = -1;
                sourceNewProc[bndI] = -1;
            }
        }
    }
}


// Subset the neighbourCell/neighbourProc fields
void Foam::fvMeshDistribute::subsetBoundaryData
(
    const fvMesh& mesh,
    const labelList& faceMap,
    const labelList& cellMap,

    const labelList& oldDistribution,
    const labelList& oldFaceOwner,
    const labelList& oldFaceNeighbour,
    const label oldInternalFaces,

    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourceNewProc,

    labelList& subFace,
    labelList& subProc,
    labelList& subNewProc
)
{
    subFace.setSize(mesh.nFaces() - mesh.nInternalFaces());
    subProc.setSize(mesh.nFaces() - mesh.nInternalFaces());
    subNewProc.setSize(mesh.nFaces() - mesh.nInternalFaces());

    forAll(subFace, newBFaceI)
    {
        label newFaceI = newBFaceI + mesh.nInternalFaces();

        label oldFaceI = faceMap[newFaceI];

        // Was oldFaceI internal face? If so which side did we get.
        if (oldFaceI < oldInternalFaces)
        {
            subFace[newBFaceI] = oldFaceI;
            subProc[newBFaceI] = Pstream::myProcNo();

            label oldOwn = oldFaceOwner[oldFaceI];
            label oldNei = oldFaceNeighbour[oldFaceI];

            if (oldOwn == cellMap[mesh.faceOwner()[newFaceI]])
            {
                // We kept the owner side. Where does the neighbour move to?
                subNewProc[newBFaceI] = oldDistribution[oldNei];
            }
            else
            {
                // We kept the neighbour side.
                subNewProc[newBFaceI] = oldDistribution[oldOwn];
            }
        }
        else
        {
            // Was boundary face. Take over boundary information
            label oldBFaceI = oldFaceI - oldInternalFaces;

            subFace[newBFaceI] = sourceFace[oldBFaceI];
            subProc[newBFaceI] = sourceProc[oldBFaceI];
            subNewProc[newBFaceI] = sourceNewProc[oldBFaceI];
        }
    }
}


// Find cells on mesh whose faceID/procID match the neighbour cell/proc of
// domainMesh. Store the matching face.
void Foam::fvMeshDistribute::findCouples
(
    const primitiveMesh& mesh,
    const labelList& sourceFace,
    const labelList& sourceProc,

    const label domain,
    const primitiveMesh& domainMesh,
    const labelList& domainFace,
    const labelList& domainProc,

    labelList& masterCoupledFaces,
    labelList& slaveCoupledFaces
)
{
    // Store domain neighbour as map so we can easily look for pair
    // with same face+proc.
    HashTable<label, labelPair, labelPair::Hash<> > map(domainFace.size());

    forAll(domainFace, bFaceI)
    {
        map.insert(labelPair(domainFace[bFaceI], domainProc[bFaceI]), bFaceI);
    }


    // Try to match mesh data.

    masterCoupledFaces.setSize(domainFace.size());
    slaveCoupledFaces.setSize(domainFace.size());
    label coupledI = 0;

    forAll(sourceFace, bFaceI)
    {
        if (sourceProc[bFaceI] != -1)
        {
            labelPair myData(sourceFace[bFaceI], sourceProc[bFaceI]);

            HashTable<label, labelPair, labelPair::Hash<> >::const_iterator
                iter = map.find(myData);

            if (iter != map.end())
            {
                label nbrBFaceI = iter();

                masterCoupledFaces[coupledI] = mesh.nInternalFaces() + bFaceI;
                slaveCoupledFaces[coupledI] =
                    domainMesh.nInternalFaces()
                  + nbrBFaceI;

                coupledI++;
            }
        }
    }

    masterCoupledFaces.setSize(coupledI);
    slaveCoupledFaces.setSize(coupledI);

    if (debug)
    {
        Pout<< "findCouples : found " << coupledI
            << " faces that will be stitched" << nl << endl;
    }
}


// Map data on boundary faces to new mesh (resulting from adding two meshes)
Foam::labelList Foam::fvMeshDistribute::mapBoundaryData
(
    const primitiveMesh& mesh,      // mesh after adding
    const mapAddedPolyMesh& map,
    const labelList& boundaryData0, // mesh before adding
    const label nInternalFaces1,
    const labelList& boundaryData1  // added mesh
)
{
    labelList newBoundaryData(mesh.nFaces() - mesh.nInternalFaces());

    forAll(boundaryData0, oldBFaceI)
    {
        label newFaceI = map.oldFaceMap()[oldBFaceI + map.nOldInternalFaces()];

        // Face still exists (is necessary?) and still boundary face
        if (newFaceI >= 0 && newFaceI >= mesh.nInternalFaces())
        {
            newBoundaryData[newFaceI - mesh.nInternalFaces()] =
                boundaryData0[oldBFaceI];
        }
    }

    forAll(boundaryData1, addedBFaceI)
    {
        label newFaceI = map.addedFaceMap()[addedBFaceI + nInternalFaces1];

        if (newFaceI >= 0 && newFaceI >= mesh.nInternalFaces())
        {
            newBoundaryData[newFaceI - mesh.nInternalFaces()] =
                boundaryData1[addedBFaceI];
        }
    }

    return newBoundaryData;
}


// Remove cells. Add all exposed faces to patch oldInternalPatchI
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::doRemoveCells
(
    const labelList& cellsToRemove,
    const label oldInternalPatchI
)
{
    // Mesh change engine
    polyTopoChange meshMod(mesh_);

    // Cell removal topo engine. Do NOT synchronize parallel since
    // we are doing a local cell removal.
    removeCells cellRemover(mesh_, false);

    // Get all exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));

    // Insert the topo changes
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        labelList(exposedFaces.size(), oldInternalPatchI),  // patch for exposed
                                                            // faces.
        meshMod
    );

    // Change the mesh. No inflation. Note: no parallel comms allowed.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, false);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }

    return map;
}


// Delete and add processor patches. Changes mesh and returns per neighbour proc
// the processor patchID.
void Foam::fvMeshDistribute::addProcPatches
(
    const labelList& neighbourNewProc,   // processor that neighbour is on
    labelList& procPatchID
)
{
    // Now use the neighbourFace/Proc to repatch the mesh. These two lists
    // contain for all current boundary faces the global patchID (for non-proc
    // patch) or the processor.

    labelList procPatchSizes(Pstream::nProcs(), 0);

    forAll(neighbourNewProc, bFaceI)
    {
        if (neighbourNewProc[bFaceI] != -1)
        {
            procPatchSizes[neighbourNewProc[bFaceI]]++;
        }
    }

    // Per neighbour processor the label of the processor patch
    procPatchID.setSize(Pstream::nProcs());

    forAll(procPatchSizes, procI)
    {
        if (procPatchSizes[procI] > 0)
        {
            const word patchName =
                "procBoundary"
              + name(Pstream::myProcNo())
              + "to"
              + name(procI);


            procPatchID[procI] = addProcPatch(patchName, procI);
            addPatchFields<volScalarField>
            (
                processorFvPatchField<scalar>::typeName
            );
            addPatchFields<volVectorField>
            (
                processorFvPatchField<vector>::typeName
            );
            addPatchFields<volSphericalTensorField>
            (
                processorFvPatchField<sphericalTensor>::typeName
            );
            addPatchFields<volSymmTensorField>
            (
                processorFvPatchField<symmTensor>::typeName
            );
            addPatchFields<volTensorField>
            (
                processorFvPatchField<tensor>::typeName
            );

            addPatchFields<surfaceScalarField>
            (
                processorFvPatchField<scalar>::typeName
            );
            addPatchFields<surfaceVectorField>
            (
                processorFvPatchField<vector>::typeName
            );
            addPatchFields<surfaceSphericalTensorField>
            (
                processorFvPatchField<sphericalTensor>::typeName
            );
            addPatchFields<surfaceSymmTensorField>
            (
                processorFvPatchField<symmTensor>::typeName
            );
            addPatchFields<surfaceTensorField>
            (
                processorFvPatchField<tensor>::typeName
            );
        }
        else
        {
            procPatchID[procI] = -1;
        }
    }
}


// Get boundary faces to be repatched. Is -1 or new patchID
Foam::labelList Foam::fvMeshDistribute::getProcBoundaryPatch
(
    const labelList& neighbourNewProc,  // new processor per boundary face
    const labelList& procPatchID        // patchID
)
{
    labelList patchIDs(neighbourNewProc);

    forAll(neighbourNewProc, bFaceI)
    {
        if (neighbourNewProc[bFaceI] != -1)
        {
            label nbrProc = neighbourNewProc[bFaceI];

            patchIDs[bFaceI] = procPatchID[nbrProc];
        }
        else
        {
            patchIDs[bFaceI] = -1;
        }
    }
    return patchIDs;
}


// Send mesh and coupling data.
void Foam::fvMeshDistribute::sendMesh
(
    const label domain,
    const fvMesh& mesh,

    const wordList& pointZoneNames,
    const wordList& faceZoneNames,
    const wordList& cellZoneNames,

    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourceNewProc,
    OSstream& toDomain
)
{
    if (debug)
    {
        Pout<< "Sending to domain " << domain << nl
            << "    nPoints:" << mesh.nPoints() << nl
            << "    nFaces:" << mesh.nFaces() << nl
            << "    nCells:" << mesh.nCells() << nl
            << "    nPatches:" << mesh.boundaryMesh().size() << nl
            << endl;
    }

    // Assume sparse, possibly overlapping point zones. Get contents
    // in merged-zone indices.
    CompactListList_dev<label> zonePoints;
    {
        const pointZoneMesh& pointZones = mesh.pointZones();

        labelList rowSizes(pointZoneNames.size(), 0);

        forAll(pointZoneNames, nameI)
        {
            label myZoneID = pointZones.findZoneID(pointZoneNames[nameI]);

            if (myZoneID != -1)
            {
                rowSizes[nameI] = pointZones[myZoneID].size();
            }
        }
        zonePoints.setSize(rowSizes);

        forAll(pointZoneNames, nameI)
        {
            label myZoneID = pointZones.findZoneID(pointZoneNames[nameI]);

            if (myZoneID != -1)
            {
                zonePoints[nameI].assign(pointZones[myZoneID]);
            }
        }
    }

    // Assume sparse, possibly overlapping face zones
    CompactListList_dev<label> zoneFaces;
    CompactListList_dev<bool> zoneFaceFlip;
    {
        const faceZoneMesh& faceZones = mesh.faceZones();

        labelList rowSizes(faceZoneNames.size(), 0);

        forAll(faceZoneNames, nameI)
        {
            label myZoneID = faceZones.findZoneID(faceZoneNames[nameI]);

            if (myZoneID != -1)
            {
                rowSizes[nameI] = faceZones[myZoneID].size();
            }
        }

        zoneFaces.setSize(rowSizes);
        zoneFaceFlip.setSize(rowSizes);

        forAll(faceZoneNames, nameI)
        {
            label myZoneID = faceZones.findZoneID(faceZoneNames[nameI]);

            if (myZoneID != -1)
            {
                zoneFaces[nameI].assign(faceZones[myZoneID]);
                zoneFaceFlip[nameI].assign(faceZones[myZoneID].flipMap());
            }
        }
    }

    // Assume sparse, possibly overlapping cell zones
    CompactListList_dev<label> zoneCells;
    {
        const cellZoneMesh& cellZones = mesh.cellZones();

        labelList rowSizes(cellZoneNames.size(), 0);

        forAll(cellZoneNames, nameI)
        {
            label myZoneID = cellZones.findZoneID(cellZoneNames[nameI]);

            if (myZoneID != -1)
            {
                rowSizes[nameI] = cellZones[myZoneID].size();
            }
        }

        zoneCells.setSize(rowSizes);

        forAll(cellZoneNames, nameI)
        {
            label myZoneID = cellZones.findZoneID(cellZoneNames[nameI]);

            if (myZoneID != -1)
            {
                zoneCells[nameI].assign(cellZones[myZoneID]);
            }
        }
    }
    ////- Assume full cell zones
    //labelList cellZoneID;
    //if (hasCellZones)
    //{
    //    cellZoneID.setSize(mesh.nCells());
    //    cellZoneID = -1;
    //
    //    const cellZoneMesh& cellZones = mesh.cellZones();
    //
    //    forAll(cellZones, zoneI)
    //    {
    //        UIndirectList<label>(cellZoneID, cellZones[zoneI]) = zoneI;
    //    }
    //}

    // Send
    toDomain
        << mesh.points()
        << CompactListList_dev<label, face>(mesh.faces())
        << mesh.faceOwner()
        << mesh.faceNeighbour()
        << mesh.boundaryMesh()

        << zonePoints
        << zoneFaces
        << zoneFaceFlip
        << zoneCells

        << sourceFace
        << sourceProc
        << sourceNewProc;


    if (debug)
    {
        Pout<< "Started sending mesh to domain " << domain
            << endl;
    }
}


// Receive mesh. Opposite of sendMesh
Foam::autoPtr<Foam::fvMesh> Foam::fvMeshDistribute::receiveMesh
(
    const label domain,
    const wordList& pointZoneNames,
    const wordList& faceZoneNames,
    const wordList& cellZoneNames,
    const Time& runTime,
    labelList& domainSourceFace,
    labelList& domainSourceProc,
    labelList& domainSourceNewProc,
    ISstream& fromNbr
)
{
    pointField domainPoints(fromNbr);
    faceList domainFaces = CompactListList_dev<label, face>(fromNbr)();
    labelList domainAllOwner(fromNbr);
    labelList domainAllNeighbour(fromNbr);
    PtrList<entry> patchEntries(fromNbr);

    CompactListList_dev<label> zonePoints(fromNbr);
    CompactListList_dev<label> zoneFaces(fromNbr);
    CompactListList_dev<bool> zoneFaceFlip(fromNbr);
    CompactListList_dev<label> zoneCells(fromNbr);

    fromNbr
        >> domainSourceFace
        >> domainSourceProc
        >> domainSourceNewProc;

    // Construct fvMesh
    autoPtr<fvMesh> domainMeshPtr
    (
        new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::NO_READ
            ),
            xferMove(domainPoints),
            xferMove(domainFaces),
            xferMove(domainAllOwner),
            xferMove(domainAllNeighbour),
            false                   // no parallel comms
        )
    );
    fvMesh& domainMesh = domainMeshPtr();

    List<polyPatch*> patches(patchEntries.size());

    forAll(patchEntries, patchI)
    {
        patches[patchI] = polyPatch::New
        (
            patchEntries[patchI].keyword(),
            patchEntries[patchI].dict(),
            patchI,
            domainMesh.boundaryMesh()
        ).ptr();
    }
    // Add patches; no parallel comms
    domainMesh.addFvPatches(patches, false);

    // Construct zones
    List<pointZone*> pZonePtrs(pointZoneNames.size());
    forAll(pZonePtrs, i)
    {
        pZonePtrs[i] = new pointZone
        (
            pointZoneNames[i],
            zonePoints[i],
            i,
            domainMesh.pointZones()
        );
    }

    List<faceZone*> fZonePtrs(faceZoneNames.size());
    forAll(fZonePtrs, i)
    {
        fZonePtrs[i] = new faceZone
        (
            faceZoneNames[i],
            zoneFaces[i],
            zoneFaceFlip[i],
            i,
            domainMesh.faceZones()
        );
    }

    List<cellZone*> cZonePtrs(cellZoneNames.size());
    forAll(cZonePtrs, i)
    {
        cZonePtrs[i] = new cellZone
        (
            cellZoneNames[i],
            zoneCells[i],
            i,
            domainMesh.cellZones()
        );
    }
    domainMesh.addZones(pZonePtrs, fZonePtrs, cZonePtrs);

    return domainMeshPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::fvMeshDistribute::fvMeshDistribute(fvMesh& mesh, const scalar mergeTol)
:
    mesh_(mesh),
    mergeTol_(mergeTol)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::fvMeshDistribute::countCells
(
    const labelList& distribution
)
{
    labelList nCells(Pstream::nProcs(), 0);
    forAll(distribution, cellI)
    {
        label newProc = distribution[cellI];

        if (newProc < 0 || newProc >= Pstream::nProcs())
        {
            FatalErrorIn("fvMeshDistribute::distribute(const labelList&)")
                << "Distribution should be in range 0.." << Pstream::nProcs()-1
                << endl
                << "At index " << cellI << " distribution:" << newProc
                << abort(FatalError);
        }
        nCells[newProc]++;
    }
    return nCells;
}


Foam::autoPtr<Foam::mapDistributePolyMesh> Foam::fvMeshDistribute::distribute
(
    const labelList& distribution
)
{
    // Some checks on distribution
    if (distribution.size() != mesh_.nCells())
    {
        FatalErrorIn("fvMeshDistribute::distribute(const labelList&)")
            << "Size of distribution:"
            << distribution.size() << " mesh nCells:" << mesh_.nCells()
            << abort(FatalError);
    }


    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Check all processors have same non-proc patches in same order.
    if (patches.checkParallelSync(true))
    {
        FatalErrorIn("fvMeshDistribute::distribute(const labelList&)")
            << "This application requires all non-processor patches"
            << " to be present in the same order on all patches" << nl
            << "followed by the processor patches (which of course are unique)."
            << nl
            << "Local patches:" << mesh_.boundaryMesh().names()
            << abort(FatalError);
    }

    // Save some data for mapping later on
    const label nOldPoints(mesh_.nPoints());
    const label nOldFaces(mesh_.nFaces());
    const label nOldCells(mesh_.nCells());
    labelList oldPatchStarts(patches.size());
    labelList oldPatchNMeshPoints(patches.size());
    forAll(patches, patchI)
    {
        oldPatchStarts[patchI] = patches[patchI].start();
        oldPatchNMeshPoints[patchI] = patches[patchI].nPoints();
    }



    // Short circuit trivial case.
    if (!Pstream::parRun())
    {
        // Collect all maps and return
        return autoPtr<mapDistributePolyMesh>
        (
            new mapDistributePolyMesh
            (
                mesh_,

                nOldPoints,
                nOldFaces,
                nOldCells,
                oldPatchStarts.xfer(),
                oldPatchNMeshPoints.xfer(),

                labelListList(1, identity(mesh_.nPoints())).xfer(),//subPointMap
                labelListList(1, identity(mesh_.nFaces())).xfer(), //subFaceMap
                labelListList(1, identity(mesh_.nCells())).xfer(), //subCellMap
                labelListList(1, identity(patches.size())).xfer(), //subPatchMap

                labelListList(1, identity(mesh_.nPoints())).xfer(),//pointMap
                labelListList(1, identity(mesh_.nFaces())).xfer(), //faceMap
                labelListList(1, identity(mesh_.nCells())).xfer(), //cellMap
                labelListList(1, identity(patches.size())).xfer()  //patchMap
            )
        );
    }


    // Collect any zone names
    const wordList pointZoneNames(mergeWordList(mesh_.pointZones().names()));
    const wordList faceZoneNames(mergeWordList(mesh_.faceZones().names()));
    const wordList cellZoneNames(mergeWordList(mesh_.cellZones().names()));



    // Local environment of all boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // A face is uniquely defined by
    //  - proc
    //  - local face no
    //
    // To glue the parts of meshes which can get sent from anywhere we
    // need to know on boundary faces what the above tuple on both sides is.
    // So we need to maintain:
    //  - original face
    //  - original processor id (= trivial)
    // For coupled boundaries (where the faces are 'duplicate') we take the
    // lowest numbered processor as the data to store.
    //
    // Additionally to create the procboundaries we need to know where the owner
    // cell on the other side now is: newNeighbourProc.
    //

    // physical boundary:
    //     sourceProc = -1
    //     sourceNewProc = -1
    //     sourceFace = patchID
    // coupled boundary:
    //     sourceProc = proc
    //     sourceNewProc = distribution[cell on proc]
    //     sourceFace = face
    labelList sourceFace;
    labelList sourceProc;
    labelList sourceNewProc;
    getNeighbourData(distribution, sourceFace, sourceProc, sourceNewProc);


    // Remove meshPhi. Since this would otherwise dissappear anyway
    // during topo changes and we have to guarantee that all the fields
    // can be sent.
    mesh_.clearOut();
    mesh_.resetMotion();

    // Get data to send. Make sure is synchronised
    const wordList volScalars(mesh_.names(volScalarField::typeName));
    checkEqualWordList("volScalarFields", volScalars);
    const wordList volVectors(mesh_.names(volVectorField::typeName));
    checkEqualWordList("volVectorFields", volVectors);
    const wordList volSphereTensors
    (
        mesh_.names(volSphericalTensorField::typeName)
    );
    checkEqualWordList("volSphericalTensorFields", volSphereTensors);
    const wordList volSymmTensors(mesh_.names(volSymmTensorField::typeName));
    checkEqualWordList("volSymmTensorFields", volSymmTensors);
    const wordList volTensors(mesh_.names(volTensorField::typeName));
    checkEqualWordList("volTensorField", volTensors);

    const wordList surfScalars(mesh_.names(surfaceScalarField::typeName));
    checkEqualWordList("surfaceScalarFields", surfScalars);
    const wordList surfVectors(mesh_.names(surfaceVectorField::typeName));
    checkEqualWordList("surfaceVectorFields", surfVectors);
    const wordList surfSphereTensors
    (
        mesh_.names(surfaceSphericalTensorField::typeName)
    );
    checkEqualWordList("surfaceSphericalTensorFields", surfSphereTensors);
    const wordList surfSymmTensors
    (
        mesh_.names(surfaceSymmTensorField::typeName)
    );
    checkEqualWordList("surfaceSymmTensorFields", surfSymmTensors);
    const wordList surfTensors(mesh_.names(surfaceTensorField::typeName));
    checkEqualWordList("surfaceTensorFields", surfTensors);




    // Find patch to temporarily put exposed and processor faces into.
    label oldInternalPatchI = findNonEmptyPatch();



    // Delete processor patches, starting from the back. Move all faces into
    // oldInternalPatchI.
    labelList repatchFaceMap;
    {
        autoPtr<mapPolyMesh> repatchMap = deleteProcPatches(oldInternalPatchI);

        // Store face map (only face ordering that changed)
        repatchFaceMap = repatchMap().faceMap();

        // Reorder all boundary face data (sourceProc, sourceFace etc.)
        labelList bFaceMap
        (
            SubList<label>
            (
                repatchMap().reverseFaceMap(),
                mesh_.nFaces() - mesh_.nInternalFaces(),
                mesh_.nInternalFaces()
            )
          - mesh_.nInternalFaces()
        );

        inplaceReorder(bFaceMap, sourceFace);
        inplaceReorder(bFaceMap, sourceProc);
        inplaceReorder(bFaceMap, sourceNewProc);
    }



    // Print a bit.
    if (debug)
    {
        Pout<< nl << "MESH WITH PROC PATCHES DELETED:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        Pout<< nl << endl;
    }



    // Maps from subsetted mesh (that is sent) back to original maps
    labelListList subCellMap(Pstream::nProcs());
    labelListList subFaceMap(Pstream::nProcs());
    labelListList subPointMap(Pstream::nProcs());
    labelListList subPatchMap(Pstream::nProcs());
    // Maps from subsetted mesh to reconstructed mesh
    labelListList constructCellMap(Pstream::nProcs());
    labelListList constructFaceMap(Pstream::nProcs());
    labelListList constructPointMap(Pstream::nProcs());
    labelListList constructPatchMap(Pstream::nProcs());




    // Find out schedule
    // ~~~~~~~~~~~~~~~~~

    labelListList nSendCells(Pstream::nProcs());
    nSendCells[Pstream::myProcNo()] = countCells(distribution);
    Pstream::gatherList(nSendCells);
    Pstream::scatterList(nSendCells);


    // Allocate buffers
    PtrList<OStringStream> sendStr(Pstream::nProcs());
    forAll(sendStr, i)
    {
        sendStr.set(i, new OStringStream(IOstream::BINARY));
    }
    //PstreamBuffers pBufs(Pstream::nonBlocking);


    // What to send to neighbouring domains
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(nSendCells[Pstream::myProcNo()], recvProc)
    {
        if
        (
            recvProc != Pstream::myProcNo()
         && nSendCells[Pstream::myProcNo()][recvProc] > 0
        )
        {
            // Send to recvProc

            if (debug)
            {
                Pout<< nl
                    << "SUBSETTING FOR DOMAIN " << recvProc
                    << " cells to send:"
                    << nSendCells[Pstream::myProcNo()][recvProc]
                    << nl << endl;
            }

            // Pstream for sending mesh and fields
            //OPstream str(Pstream::blocking, recvProc);
            //UOPstream str(recvProc, pBufs);

            // Mesh subsetting engine
            fvMeshSubset subsetter(mesh_);

            // Subset the cells of the current domain.
            subsetter.setLargeCellSubset
            (
                distribution,
                recvProc,
                oldInternalPatchI,  // oldInternalFaces patch
                false               // no parallel sync
            );

            subCellMap[recvProc] = subsetter.cellMap();
            subFaceMap[recvProc] = renumber
            (
                repatchFaceMap,
                subsetter.faceMap()
            );
            subPointMap[recvProc] = subsetter.pointMap();
            subPatchMap[recvProc] = subsetter.patchMap();


            // Subset the boundary fields (owner/neighbour/processor)
            labelList procSourceFace;
            labelList procSourceProc;
            labelList procSourceNewProc;

            subsetBoundaryData
            (
                subsetter.subMesh(),
                subsetter.faceMap(),        // from subMesh to mesh
                subsetter.cellMap(),        //      ,,      ,,

                distribution,               // old mesh distribution
                mesh_.faceOwner(),          // old owner
                mesh_.faceNeighbour(),
                mesh_.nInternalFaces(),

                sourceFace,
                sourceProc,
                sourceNewProc,

                procSourceFace,
                procSourceProc,
                procSourceNewProc
            );



            // Send to neighbour
            sendMesh
            (
                recvProc,
                subsetter.subMesh(),

                pointZoneNames,
                faceZoneNames,
                cellZoneNames,

                procSourceFace,
                procSourceProc,
                procSourceNewProc,
                sendStr[recvProc]
            );
            sendFields<volScalarField>
            (
                recvProc,
                volScalars,
                subsetter,
                sendStr[recvProc]
            );
            sendFields<volVectorField>
            (
                recvProc,
                volVectors,
                subsetter,
                sendStr[recvProc]
            );
            sendFields<volSphericalTensorField>
            (
                recvProc,
                volSphereTensors,
                subsetter,
                sendStr[recvProc]
            );
            sendFields<volSymmTensorField>
            (
                recvProc,
                volSymmTensors,
                subsetter,
                sendStr[recvProc]
            );
            sendFields<volTensorField>
            (
                recvProc,
                volTensors,
                subsetter,
                sendStr[recvProc]
            );

            sendFields<surfaceScalarField>
            (
                recvProc,
                surfScalars,
                subsetter,
                sendStr[recvProc]
            );
            sendFields<surfaceVectorField>
            (
                recvProc,
                surfVectors,
                subsetter,
                sendStr[recvProc]
            );
            sendFields<surfaceSphericalTensorField>
            (
                recvProc,
                surfSphereTensors,
                subsetter,
                sendStr[recvProc]
            );
            sendFields<surfaceSymmTensorField>
            (
                recvProc,
                surfSymmTensors,
                subsetter,
                sendStr[recvProc]
            );
            sendFields<surfaceTensorField>
            (
                recvProc,
                surfTensors,
                subsetter,
                sendStr[recvProc]
            );
        }
    }


    // Start sending&receiving from buffers
    //pBufs.finishedSends();

    // get the data.
    PtrList<IStringStream> recvStr(Pstream::nProcs());
    {
        List<List<char> > sendBufs(sendStr.size());
        forAll(sendStr, procI)
        {
            string contents = sendStr[procI].str();
            const char* ptr = contents.data();

            sendBufs[procI].setSize(contents.size());
            forAll(sendBufs[procI], i)
            {
                sendBufs[procI][i] = *ptr++;
            }
            // Clear OStringStream
            sendStr.set(procI, NULL);
        }

        // Transfer sendBufs into recvBufs
        List<List<char> > recvBufs(Pstream::nProcs());
        labelListList sizes(Pstream::nProcs());
        exchange<List<char>, char>(sendBufs, recvBufs, sizes);

        forAll(recvStr, procI)
        {
            string contents(recvBufs[procI].begin(), recvBufs[procI].size());
            recvStr.set
            (
                procI,
                new IStringStream(contents, IOstream::BINARY)
            );
        }
    }


    // Subset the part that stays
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Save old mesh maps before changing mesh
        const labelList oldFaceOwner(mesh_.faceOwner());
        const labelList oldFaceNeighbour(mesh_.faceNeighbour());
        const label oldInternalFaces = mesh_.nInternalFaces();

        // Remove cells.
        autoPtr<mapPolyMesh> subMap
        (
            doRemoveCells
            (
                select(false, distribution, Pstream::myProcNo()),
                oldInternalPatchI
            )
        );

        // Addressing from subsetted mesh
        subCellMap[Pstream::myProcNo()] = subMap().cellMap();
        subFaceMap[Pstream::myProcNo()] = renumber
        (
            repatchFaceMap,
            subMap().faceMap()
        );
        subPointMap[Pstream::myProcNo()] = subMap().pointMap();
        subPatchMap[Pstream::myProcNo()] = identity(patches.size());

        // Initialize all addressing into current mesh
        constructCellMap[Pstream::myProcNo()] = identity(mesh_.nCells());
        constructFaceMap[Pstream::myProcNo()] = identity(mesh_.nFaces());
        constructPointMap[Pstream::myProcNo()] = identity(mesh_.nPoints());
        constructPatchMap[Pstream::myProcNo()] = identity(patches.size());

        // Subset the mesh data: neighbourCell/neighbourProc
        // fields
        labelList domainSourceFace;
        labelList domainSourceProc;
        labelList domainSourceNewProc;

        subsetBoundaryData
        (
            mesh_,                          // new mesh
            subMap().faceMap(),             // from new to original mesh
            subMap().cellMap(),

            distribution,                   // distribution before subsetting
            oldFaceOwner,                   // owner before subsetting
            oldFaceNeighbour,               // neighbour        ,,
            oldInternalFaces,               // nInternalFaces   ,,

            sourceFace,
            sourceProc,
            sourceNewProc,

            domainSourceFace,
            domainSourceProc,
            domainSourceNewProc
        );

        sourceFace.transfer(domainSourceFace);
        sourceProc.transfer(domainSourceProc);
        sourceNewProc.transfer(domainSourceNewProc);
    }


    // Print a bit.
    if (debug)
    {
        Pout<< nl << "STARTING MESH:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        Pout<< nl << endl;
    }



    // Receive and add what was sent
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(nSendCells, sendProc)
    {
        // Did processor sendProc send anything to me?
        if
        (
            sendProc != Pstream::myProcNo()
         && nSendCells[sendProc][Pstream::myProcNo()] > 0
        )
        {
            if (debug)
            {
                Pout<< nl
                    << "RECEIVING FROM DOMAIN " << sendProc
                    << " cells to receive:"
                    << nSendCells[sendProc][Pstream::myProcNo()]
                    << nl << endl;
            }


            // Pstream for receiving mesh and fields
            //UIPstream str(sendProc, pBufs);


            // Receive from sendProc
            labelList domainSourceFace;
            labelList domainSourceProc;
            labelList domainSourceNewProc;

            autoPtr<fvMesh> domainMeshPtr;
            PtrList<volScalarField> vsf;
            PtrList<volVectorField> vvf;
            PtrList<volSphericalTensorField> vsptf;
            PtrList<volSymmTensorField> vsytf;
            PtrList<volTensorField> vtf;
            PtrList<surfaceScalarField> ssf;
            PtrList<surfaceVectorField> svf;
            PtrList<surfaceSphericalTensorField> ssptf;
            PtrList<surfaceSymmTensorField> ssytf;
            PtrList<surfaceTensorField> stf;

            // Opposite of sendMesh
            {
                domainMeshPtr = receiveMesh
                (
                    sendProc,
                    pointZoneNames,
                    faceZoneNames,
                    cellZoneNames,

                    const_cast<Time&>(mesh_.time()),
                    domainSourceFace,
                    domainSourceProc,
                    domainSourceNewProc,
                    recvStr[sendProc]
                );
                fvMesh& domainMesh = domainMeshPtr();

                // Receive fields. Read as single dictionary because
                // of problems reading consecutive fields from single stream.
                dictionary fieldDicts(recvStr[sendProc]);

                receiveFields<volScalarField>
                (
                    sendProc,
                    volScalars,
                    domainMesh,
                    vsf,
                    fieldDicts.subDict(volScalarField::typeName)
                );
                receiveFields<volVectorField>
                (
                    sendProc,
                    volVectors,
                    domainMesh,
                    vvf,
                    fieldDicts.subDict(volVectorField::typeName)
                );
                receiveFields<volSphericalTensorField>
                (
                    sendProc,
                    volSphereTensors,
                    domainMesh,
                    vsptf,
                    fieldDicts.subDict(volSphericalTensorField::typeName)
                );
                receiveFields<volSymmTensorField>
                (
                    sendProc,
                    volSymmTensors,
                    domainMesh,
                    vsytf,
                    fieldDicts.subDict(volSymmTensorField::typeName)
                );
                receiveFields<volTensorField>
                (
                    sendProc,
                    volTensors,
                    domainMesh,
                    vtf,
                    fieldDicts.subDict(volTensorField::typeName)
                );

                receiveFields<surfaceScalarField>
                (
                    sendProc,
                    surfScalars,
                    domainMesh,
                    ssf,
                    fieldDicts.subDict(surfaceScalarField::typeName)
                );
                receiveFields<surfaceVectorField>
                (
                    sendProc,
                    surfVectors,
                    domainMesh,
                    svf,
                    fieldDicts.subDict(surfaceVectorField::typeName)
                );
                receiveFields<surfaceSphericalTensorField>
                (
                    sendProc,
                    surfSphereTensors,
                    domainMesh,
                    ssptf,
                    fieldDicts.subDict(surfaceSphericalTensorField::typeName)
                );
                receiveFields<surfaceSymmTensorField>
                (
                    sendProc,
                    surfSymmTensors,
                    domainMesh,
                    ssytf,
                    fieldDicts.subDict(surfaceSymmTensorField::typeName)
                );
                receiveFields<surfaceTensorField>
                (
                    sendProc,
                    surfTensors,
                    domainMesh,
                    stf,
                    fieldDicts.subDict(surfaceTensorField::typeName)
                );
            }
            const fvMesh& domainMesh = domainMeshPtr();


            constructCellMap[sendProc] = identity(domainMesh.nCells());
            constructFaceMap[sendProc] = identity(domainMesh.nFaces());
            constructPointMap[sendProc] = identity(domainMesh.nPoints());
            constructPatchMap[sendProc] =
                identity(domainMesh.boundaryMesh().size());


            // Print a bit.
            if (debug)
            {
                Pout<< nl << "RECEIVED MESH FROM:" << sendProc << endl;
                printMeshInfo(domainMesh);
                printFieldInfo<volScalarField>(domainMesh);
                printFieldInfo<volVectorField>(domainMesh);
                printFieldInfo<volSphericalTensorField>(domainMesh);
                printFieldInfo<volSymmTensorField>(domainMesh);
                printFieldInfo<volTensorField>(domainMesh);
                printFieldInfo<surfaceScalarField>(domainMesh);
                printFieldInfo<surfaceVectorField>(domainMesh);
                printFieldInfo<surfaceSphericalTensorField>(domainMesh);
                printFieldInfo<surfaceSymmTensorField>(domainMesh);
                printFieldInfo<surfaceTensorField>(domainMesh);
            }


            // Now this mesh we received (from sendProc) needs to be merged
            // with the current mesh. On the current mesh we have for all
            // boundaryfaces the original face and processor. See if we can
            // match these up to the received domainSourceFace and
            // domainSourceProc.
            labelList masterCoupledFaces;
            labelList slaveCoupledFaces;
            findCouples
            (
                mesh_,

                sourceFace,
                sourceProc,

                sendProc,
                domainMesh,
                domainSourceFace,
                domainSourceProc,

                masterCoupledFaces,
                slaveCoupledFaces
            );

            // Generate additional coupling info (points, edges) from
            // faces-that-match
            faceCoupleInfo couples
            (
                mesh_,
                masterCoupledFaces,
                domainMesh,
                slaveCoupledFaces,
                mergeTol_,              // merge tolerance
                true,                   // faces align
                true,                   // couples are ordered already
                false
            );


            // Add domainMesh to mesh
            // ~~~~~~~~~~~~~~~~~~~~~~

            autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
            (
                mesh_,
                domainMesh,
                couples,
                false           // no parallel comms
            );

            // Update mesh data: sourceFace,sourceProc for added
            // mesh.

            sourceFace =
                mapBoundaryData
                (
                    mesh_,
                    map(),
                    sourceFace,
                    domainMesh.nInternalFaces(),
                    domainSourceFace
                );
            sourceProc =
                mapBoundaryData
                (
                    mesh_,
                    map(),
                    sourceProc,
                    domainMesh.nInternalFaces(),
                    domainSourceProc
                );
            sourceNewProc =
                mapBoundaryData
                (
                    mesh_,
                    map(),
                    sourceNewProc,
                    domainMesh.nInternalFaces(),
                    domainSourceNewProc
                );

            // Update all addressing so xxProcAddressing points to correct item
            // in masterMesh.
            const labelList& oldCellMap = map().oldCellMap();
            const labelList& oldFaceMap = map().oldFaceMap();
            const labelList& oldPointMap = map().oldPointMap();
            const labelList& oldPatchMap = map().oldPatchMap();

            forAll(constructPatchMap, procI)
            {
                if (procI != sendProc && constructPatchMap[procI].size())
                {
                    // Processor already in mesh (either myProcNo or received)
                    inplaceRenumber(oldCellMap, constructCellMap[procI]);
                    inplaceRenumber(oldFaceMap, constructFaceMap[procI]);
                    inplaceRenumber(oldPointMap, constructPointMap[procI]);
                    inplaceRenumber(oldPatchMap, constructPatchMap[procI]);
                }
            }

            // Added processor
            inplaceRenumber(map().addedCellMap(), constructCellMap[sendProc]);
            inplaceRenumber(map().addedFaceMap(), constructFaceMap[sendProc]);
            inplaceRenumber(map().addedPointMap(), constructPointMap[sendProc]);
            inplaceRenumber(map().addedPatchMap(), constructPatchMap[sendProc]);

            if (debug)
            {
                Pout<< nl << "MERGED MESH FROM:" << sendProc << endl;
                printMeshInfo(mesh_);
                printFieldInfo<volScalarField>(mesh_);
                printFieldInfo<volVectorField>(mesh_);
                printFieldInfo<volSphericalTensorField>(mesh_);
                printFieldInfo<volSymmTensorField>(mesh_);
                printFieldInfo<volTensorField>(mesh_);
                printFieldInfo<surfaceScalarField>(mesh_);
                printFieldInfo<surfaceVectorField>(mesh_);
                printFieldInfo<surfaceSphericalTensorField>(mesh_);
                printFieldInfo<surfaceSymmTensorField>(mesh_);
                printFieldInfo<surfaceTensorField>(mesh_);
                Pout<< nl << endl;
            }
        }
    }


    // Print a bit.
    if (debug)
    {
        Pout<< nl << "REDISTRIBUTED MESH:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        Pout<< nl << endl;
    }



    // Add processorPatches
    // ~~~~~~~~~~~~~~~~~~~~

    // Per neighbour processor the patchID to it (or -1).
    labelList procPatchID;

    // Add processor patches.
    addProcPatches(sourceNewProc, procPatchID);

    // Put faces into correct patch. Note that we now have proper
    // processorPolyPatches again so repatching will take care of coupled face
    // ordering.

    // Get boundary faces to be repatched. Is -1 or new patchID
    labelList newPatchID
    (
        getProcBoundaryPatch
        (
            sourceNewProc,
            procPatchID
        )
    );

    // Change patches. Since this might change ordering of coupled faces
    // we also need to adapt our constructMaps.
    repatch(newPatchID, constructFaceMap);

    // See if any geometrically shared points need to be merged. Note: does
    // parallel comms.
    mergeSharedPoints(constructPointMap);

    // Bit of hack: processorFvPatchField does not get reset since created
    // from nothing so explicitly reset.
    initPatchFields<volScalarField, processorFvPatchField<scalar> >
    (
        pTraits<scalar>::zero
    );
    initPatchFields<volVectorField, processorFvPatchField<vector> >
    (
        pTraits<vector>::zero
    );
    initPatchFields
    <
        volSphericalTensorField,
        processorFvPatchField<sphericalTensor>
    >
    (
        pTraits<sphericalTensor>::zero
    );
    initPatchFields<volSymmTensorField, processorFvPatchField<symmTensor> >
    (
        pTraits<symmTensor>::zero
    );
    initPatchFields<volTensorField, processorFvPatchField<tensor> >
    (
        pTraits<tensor>::zero
    );
    initPatchFields<surfaceScalarField, processorFvsPatchField<scalar> >
    (
        pTraits<scalar>::zero
    );
    initPatchFields<surfaceVectorField, processorFvsPatchField<vector> >
    (
        pTraits<vector>::zero
    );
    initPatchFields
    <
        surfaceSphericalTensorField,
        processorFvsPatchField<sphericalTensor>
    >
    (
        pTraits<sphericalTensor>::zero
    );
    initPatchFields
    <
        surfaceSymmTensorField,
        processorFvsPatchField<symmTensor>
    >
    (
        pTraits<symmTensor>::zero
    );
    initPatchFields<surfaceTensorField, processorFvsPatchField<tensor> >
    (
        pTraits<tensor>::zero
    );


    mesh_.setInstance(mesh_.time().timeName());


    // Print a bit
    if (debug)
    {
        Pout<< nl << "FINAL MESH:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        Pout<< nl << endl;
    }

    // Collect all maps and return
    return autoPtr<mapDistributePolyMesh>
    (
        new mapDistributePolyMesh
        (
            mesh_,

            nOldPoints,
            nOldFaces,
            nOldCells,
            oldPatchStarts,
            oldPatchNMeshPoints,

            subPointMap,
            subFaceMap,
            subCellMap,
            subPatchMap,

            constructPointMap,
            constructFaceMap,
            constructCellMap,
            constructPatchMap,
            true                // reuse storage
        )
    );
}


// ************************************************************************* //
