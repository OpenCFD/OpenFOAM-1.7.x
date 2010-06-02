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

\*---------------------------------------------------------------------------*/

#include "FaceCellWave.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "OPstream.H"
#include "IPstream.H"
#include "PstreamReduceOps.H"
#include "debug.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template <class Type>
const Foam::scalar Foam::FaceCellWave<Type>::geomTol_ = 1e-6;

template <class Type>
Foam::scalar Foam::FaceCellWave<Type>::propagationTol_ = 0.01;

// Write to ostream
template <class Type>
Foam::Ostream& Foam::FaceCellWave<Type>::writeFaces
(
    const label nFaces,
    const labelList& faceLabels,
    const List<Type>& faceInfo,
    Ostream& os
)
{
    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII)
    {
        os << nFaces;

        for(label i = 0; i < nFaces; i++)
        {
            os << ' ' << faceLabels[i];
        }
        for(label i = 0; i < nFaces; i++)
        {
            os << ' ' << faceInfo[i];
        }
    }
    else
    {
        os << nFaces;

        for(label i = 0; i < nFaces; i++)
        {
            os << faceLabels[i];
        }
        for(label i = 0; i < nFaces; i++)
        {
            os << faceInfo[i];
        }
    }
    return os;
}


// Read from istream
template <class Type>
Foam::Istream& Foam::FaceCellWave<Type>::readFaces
(
    label& nFaces,
    labelList& faceLabels,
    List<Type>& faceInfo,
    Istream& is
)
{
    is >> nFaces;

    for(label i = 0; i < nFaces; i++)
    {
        is >> faceLabels[i];
    }
    for(label i = 0; i < nFaces; i++)
    {
        is >> faceInfo[i];
    }
    return is;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update info for cellI, at position pt, with information from
// neighbouring face/cell.
// Updates:
//      - changedCell_, changedCells_, nChangedCells_,
//      - statistics: nEvals_, nUnvisitedCells_
template <class Type>
bool Foam::FaceCellWave<Type>::updateCell
(
    const label cellI,
    const label neighbourFaceI,
    const Type& neighbourInfo,
    const scalar tol,
    Type& cellInfo
)
{
    nEvals_++;

    bool wasValid = cellInfo.valid();

    bool propagate =
        cellInfo.updateCell
        (
            mesh_,
            cellI,
            neighbourFaceI,
            neighbourInfo,
            tol
        );

    if (propagate)
    {
        if (!changedCell_[cellI])
        {
            changedCell_[cellI] = true;
            changedCells_[nChangedCells_++] = cellI;
        }
    }

    if (!wasValid && cellInfo.valid())
    {
        --nUnvisitedCells_;
    }

    return propagate;
}


// Update info for faceI, at position pt, with information from
// neighbouring face/cell.
// Updates:
//      - changedFace_, changedFaces_, nChangedFaces_,
//      - statistics: nEvals_, nUnvisitedFaces_
template <class Type>
bool Foam::FaceCellWave<Type>::updateFace
(
    const label faceI,
    const label neighbourCellI,
    const Type& neighbourInfo,
    const scalar tol,
    Type& faceInfo
)
{
    nEvals_++;

    bool wasValid = faceInfo.valid();

    bool propagate =
        faceInfo.updateFace
        (
            mesh_,
            faceI,
            neighbourCellI,
            neighbourInfo,
            tol
        );

    if (propagate)
    {
        if (!changedFace_[faceI])
        {
            changedFace_[faceI] = true;
            changedFaces_[nChangedFaces_++] = faceI;
        }
    }

    if (!wasValid && faceInfo.valid())
    {
        --nUnvisitedFaces_;
    }

    return propagate;
}


// Update info for faceI, at position pt, with information from
// same face.
// Updates:
//      - changedFace_, changedFaces_, nChangedFaces_,
//      - statistics: nEvals_, nUnvisitedFaces_
template <class Type>
bool Foam::FaceCellWave<Type>::updateFace
(
    const label faceI,
    const Type& neighbourInfo,
    const scalar tol,
    Type& faceInfo
)
{
    nEvals_++;

    bool wasValid = faceInfo.valid();

    bool propagate =
        faceInfo.updateFace
        (
            mesh_,
            faceI,
            neighbourInfo,
            tol
        );

    if (propagate)
    {
        if (!changedFace_[faceI])
        {
            changedFace_[faceI] = true;
            changedFaces_[nChangedFaces_++] = faceI;
        }
    }

    if (!wasValid && faceInfo.valid())
    {
        --nUnvisitedFaces_;
    }

    return propagate;
}


// For debugging: check status on both sides of cyclic
template <class Type>
void Foam::FaceCellWave<Type>::checkCyclic(const polyPatch& patch) const
{
    label cycOffset = patch.size()/2;

    for(label patchFaceI = 0; patchFaceI < cycOffset; patchFaceI++)
    {
        label i1 = patch.start() + patchFaceI;
        label i2 = i1 + cycOffset;

        if (!allFaceInfo_[i1].sameGeometry(mesh_, allFaceInfo_[i2], geomTol_))
        {
            FatalErrorIn("FaceCellWave<Type>::checkCyclic(const polyPatch&)")
                << "problem: i:" << i1 << "  otheri:" << i2
                << "   faceInfo:" << allFaceInfo_[i1]
                << "   otherfaceInfo:" << allFaceInfo_[i2]
                << abort(FatalError);
        }

        if (changedFace_[i1] != changedFace_[i2])
        {
            FatalErrorIn("FaceCellWave<Type>::checkCyclic(const polyPatch&)")
                << " problem: i:" << i1 << "  otheri:" << i2
                << "   faceInfo:" << allFaceInfo_[i1]
                << "   otherfaceInfo:" << allFaceInfo_[i2]
                << "   changedFace:" << changedFace_[i1]
                << "   otherchangedFace:" << changedFace_[i2]
                << abort(FatalError);
        }
    }
}


// Check if has cyclic patches
template <class Type>
bool Foam::FaceCellWave<Type>::hasCyclicPatch() const
{
    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (isA<cyclicPolyPatch>(mesh_.boundaryMesh()[patchI]))
        {
            return true;
        }
    }
    return false;
}


// Copy face information into member data
template <class Type>
void Foam::FaceCellWave<Type>::setFaceInfo
(
    const labelList& changedFaces,
    const List<Type>& changedFacesInfo
)
{
    forAll(changedFaces, changedFaceI)
    {
        label faceI = changedFaces[changedFaceI];

        bool wasValid = allFaceInfo_[faceI].valid();

        // Copy info for faceI
        allFaceInfo_[faceI] = changedFacesInfo[changedFaceI];

        // Maintain count of unset faces
        if (!wasValid && allFaceInfo_[faceI].valid())
        {
            --nUnvisitedFaces_;
        }

        // Mark faceI as changed, both on list and on face itself.

        changedFace_[faceI] = true;
        changedFaces_[nChangedFaces_++] = faceI;
    }
}


// Merge face information into member data
template <class Type>
void Foam::FaceCellWave<Type>::mergeFaceInfo
(
    const polyPatch& patch,
    const label nFaces,
    const labelList& changedFaces,
    const List<Type>& changedFacesInfo,
    const bool
)
{
    for(label changedFaceI = 0; changedFaceI < nFaces; changedFaceI++)
    {
        const Type& neighbourWallInfo = changedFacesInfo[changedFaceI];
        label patchFaceI = changedFaces[changedFaceI];

        label meshFaceI = patch.start() + patchFaceI;

        Type& currentWallInfo = allFaceInfo_[meshFaceI];

        if (currentWallInfo != neighbourWallInfo)
        {
            updateFace
            (
                meshFaceI,
                neighbourWallInfo,
                propagationTol_,
                currentWallInfo
            );
        }
    }
}


// Construct compact patchFace change arrays for a (slice of a) single patch.
// changedPatchFaces in local patch numbering.
// Return length of arrays.
template <class Type>
Foam::label Foam::FaceCellWave<Type>::getChangedPatchFaces
(
    const polyPatch& patch,
    const label startFaceI,
    const label nFaces,
    labelList& changedPatchFaces,
    List<Type>& changedPatchFacesInfo
) const
{
    label nChangedPatchFaces = 0;

    for(label i = 0; i < nFaces; i++)
    {
        label patchFaceI = i + startFaceI;

        label meshFaceI = patch.start() + patchFaceI;

        if (changedFace_[meshFaceI])
        {
            changedPatchFaces[nChangedPatchFaces] = patchFaceI;
            changedPatchFacesInfo[nChangedPatchFaces] = allFaceInfo_[meshFaceI];
            nChangedPatchFaces++;
        }
    }
    return nChangedPatchFaces;
}


// Handle leaving domain. Implementation referred to Type
template <class Type>
void Foam::FaceCellWave<Type>::leaveDomain
(
    const polyPatch& patch,
    const label nFaces,
    const labelList& faceLabels,
    List<Type>& faceInfo
) const
{
    const vectorField& fc = mesh_.faceCentres();

    for(label i = 0; i < nFaces; i++)
    {
        label patchFaceI = faceLabels[i];

        label meshFaceI = patch.start() + patchFaceI;
        faceInfo[i].leaveDomain(mesh_, patch, patchFaceI, fc[meshFaceI]);
    }
}


// Handle entering domain. Implementation referred to Type
template <class Type>
void Foam::FaceCellWave<Type>::enterDomain
(
    const polyPatch& patch,
    const label nFaces,
    const labelList& faceLabels,
    List<Type>& faceInfo
) const
{
    const vectorField& fc = mesh_.faceCentres();

    for(label i = 0; i < nFaces; i++)
    {
        label patchFaceI = faceLabels[i];

        label meshFaceI = patch.start() + patchFaceI;
        faceInfo[i].enterDomain(mesh_, patch, patchFaceI, fc[meshFaceI]);
    }
}


// Transform. Implementation referred to Type
template <class Type>
void Foam::FaceCellWave<Type>::transform
(
    const tensorField& rotTensor,
    const label nFaces,
    List<Type>& faceInfo
)
{
    if (rotTensor.size() == 1)
    {
        const tensor& T = rotTensor[0];

        for(label faceI = 0; faceI < nFaces; faceI++)
        {
            faceInfo[faceI].transform(mesh_, T);
        }
    }
    else
    {
        for(label faceI = 0; faceI < nFaces; faceI++)
        {
            faceInfo[faceI].transform(mesh_, rotTensor[faceI]);
        }
    }
}


// Send face info to neighbour.
template <class Type>
void Foam::FaceCellWave<Type>::sendPatchInfo
(
    const label neighbour,
    const label nFaces,
    const labelList& faceLabels,
    const List<Type>& faceInfo
) const
{
    OPstream toNeighbour(Pstream::blocking, neighbour);

    writeFaces(nFaces, faceLabels, faceInfo, toNeighbour);
}


// Receive face info from neighbour
template <class Type>
Foam::label Foam::FaceCellWave<Type>::receivePatchInfo
(
    const label neighbour,
    labelList& faceLabels,
    List<Type>& faceInfo
) const
{
    IPstream fromNeighbour(Pstream::blocking, neighbour);

    label nFaces = 0;
    readFaces(nFaces, faceLabels, faceInfo, fromNeighbour);

    return nFaces;
}


// Offset mesh face. Used for transferring from one cyclic half to the other.
template <class Type>
void Foam::FaceCellWave<Type>::offset
(
    const polyPatch&,
    const label cycOffset,
    const label nFaces,
    labelList& faces
)
{
    for(label faceI = 0; faceI < nFaces; faceI++)
    {
        faces[faceI] += cycOffset;
    }
}


// Tranfer all the information to/from neighbouring processors
template <class Type>
void Foam::FaceCellWave<Type>::handleProcPatches()
{
    // Send all

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(patch))
        {
            // Allocate buffers
            label nSendFaces;
            labelList sendFaces(patch.size());
            List<Type> sendFacesInfo(patch.size());

            // Determine which faces changed on current patch
            nSendFaces = getChangedPatchFaces
            (
                patch,
                0,
                patch.size(),
                sendFaces,
                sendFacesInfo
            );

            // Adapt wallInfo for leaving domain
            leaveDomain
            (
                patch,
                nSendFaces,
                sendFaces,
                sendFacesInfo
            );

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patch);

            if (debug)
            {
                Pout<< " Processor patch " << patchI << ' ' << patch.name()
                    << " communicating with " << procPatch.neighbProcNo()
                    << "  Sending:" << nSendFaces
                    << endl;
            }

            sendPatchInfo
            (
                procPatch.neighbProcNo(),
                nSendFaces,
                sendFaces,
                sendFacesInfo
            );
        }
    }

    // Receive all

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(patch))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patch);

            // Allocate buffers
            label nReceiveFaces;
            labelList receiveFaces(patch.size());
            List<Type> receiveFacesInfo(patch.size());

            nReceiveFaces = receivePatchInfo
            (
                procPatch.neighbProcNo(),
                receiveFaces,
                receiveFacesInfo
            );

            if (debug)
            {
                Pout<< " Processor patch " << patchI << ' ' << patch.name()
                    << " communicating with " << procPatch.neighbProcNo()
                    << "  Receiving:" << nReceiveFaces
                    << endl;
            }

            // Apply transform to received data for non-parallel planes
            if (!procPatch.parallel())
            {
                transform
                (
                    procPatch.reverseT(),
                    nReceiveFaces,
                    receiveFacesInfo
                );
            }

            // Adapt wallInfo for entering domain
            enterDomain
            (
                patch,
                nReceiveFaces,
                receiveFaces,
                receiveFacesInfo
            );

            // Merge received info
            mergeFaceInfo
            (
                patch,
                nReceiveFaces,
                receiveFaces,
                receiveFacesInfo,
                procPatch.parallel()
            );
        }
    }
}


// Transfer information across cyclic halves.
// Note: Cyclic is two patches in one: one side from 0..size/2 and other
// side from size/2 .. size.
template <class Type>
void Foam::FaceCellWave<Type>::handleCyclicPatches()
{
    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];

        if (isA<cyclicPolyPatch>(patch))
        {
            label halfSize = patch.size()/2;

            // Allocate buffers
            label nSendFaces;
            labelList sendFaces(halfSize);
            List<Type> sendFacesInfo(halfSize);

            label nReceiveFaces;
            labelList receiveFaces(halfSize);
            List<Type> receiveFacesInfo(halfSize);

            // Half1: Determine which faces changed. Use sendFaces for storage
            nSendFaces = getChangedPatchFaces
            (
                patch,
                0,
                halfSize,
                sendFaces,
                sendFacesInfo
            );

            // Half2: Determine which faces changed. Use receiveFaces_  ,,
            nReceiveFaces = getChangedPatchFaces
            (
                patch,
                halfSize,
                halfSize,
                receiveFaces,
                receiveFacesInfo
            );

            //Info<< "Half1:" << endl;
            //writeFaces(nSendFaces, sendFaces, sendFacesInfo, Info);
            //Info<< endl;
            //
            //Info<< "Half2:" << endl;
            //writeFaces(nReceiveFaces, receiveFaces, receiveFacesInfo, Info);
            //Info<< endl;


            // Half1: Adapt wallInfo for leaving domain
            leaveDomain
            (
                patch,
                nSendFaces,
                sendFaces,
                sendFacesInfo
            );
            // Half2: Adapt wallInfo for leaving domain
            leaveDomain
            (
                patch,
                nReceiveFaces,
                receiveFaces,
                receiveFacesInfo
            );

            // Half1: 'transfer' to other side by offsetting patchFaceI
            offset(patch, halfSize, nSendFaces, sendFaces);

            // Half2: 'transfer' to other side
            offset(patch, -halfSize, nReceiveFaces, receiveFaces);

            // Apply rotation for non-parallel planes
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patch);

            if (!cycPatch.parallel())
            {
                // sendFaces = received data from half1
                transform
                (
                    cycPatch.forwardT(),
                    nSendFaces,
                    sendFacesInfo
                );

                // receiveFaces = received data from half2
                transform
                (
                    cycPatch.reverseT(),
                    nReceiveFaces,
                    receiveFacesInfo
                );
            }

            if (debug)
            {
                Pout<< " Cyclic patch " << patchI << ' ' << patch.name()
                    << "  Changed on first half : " << nSendFaces
                    << "  Changed on second half : " << nReceiveFaces
                    << endl;
            }

            // Half1: Adapt wallInfo for entering domain
            enterDomain
            (
                patch,
                nSendFaces,
                sendFaces,
                sendFacesInfo
            );

            // Half2: Adapt wallInfo for entering domain
            enterDomain
            (
                patch,
                nReceiveFaces,
                receiveFaces,
                receiveFacesInfo
            );

            // Merge into global storage
            mergeFaceInfo
            (
                patch,
                nSendFaces,
                sendFaces,
                sendFacesInfo,
                cycPatch.parallel()
            );
            // Merge into global storage
            mergeFaceInfo
            (
                patch,
                nReceiveFaces,
                receiveFaces,
                receiveFacesInfo,
                cycPatch.parallel()
            );

            if (debug)
            {
                checkCyclic(patch);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Set up only. Use setFaceInfo and iterate() to do actual calculation.
template <class Type>
Foam::FaceCellWave<Type>::FaceCellWave
(
    const polyMesh& mesh,
    UList<Type>& allFaceInfo,
    UList<Type>& allCellInfo
)
:
    mesh_(mesh),
    allFaceInfo_(allFaceInfo),
    allCellInfo_(allCellInfo),
    changedFace_(mesh_.nFaces(), false),
    changedFaces_(mesh_.nFaces()),
    nChangedFaces_(0),
    changedCell_(mesh_.nCells(), false),
    changedCells_(mesh_.nCells()),
    nChangedCells_(0),
    hasCyclicPatches_(hasCyclicPatch()),
    nEvals_(0),
    nUnvisitedCells_(mesh_.nCells()),
    nUnvisitedFaces_(mesh_.nFaces()),
    iter_(0)
{}


// Iterate, propagating changedFacesInfo across mesh, until no change (or
// maxIter reached). Initial cell values specified.
template <class Type>
Foam::FaceCellWave<Type>::FaceCellWave
(
    const polyMesh& mesh,
    const labelList& changedFaces,
    const List<Type>& changedFacesInfo,
    UList<Type>& allFaceInfo,
    UList<Type>& allCellInfo,
    const label maxIter
)
:
    mesh_(mesh),
    allFaceInfo_(allFaceInfo),
    allCellInfo_(allCellInfo),
    changedFace_(mesh_.nFaces(), false),
    changedFaces_(mesh_.nFaces()),
    nChangedFaces_(0),
    changedCell_(mesh_.nCells(), false),
    changedCells_(mesh_.nCells()),
    nChangedCells_(0),
    hasCyclicPatches_(hasCyclicPatch()),
    nEvals_(0),
    nUnvisitedCells_(mesh_.nCells()),
    nUnvisitedFaces_(mesh_.nFaces()),
    iter_(0)
{
    // Copy initial changed faces data
    setFaceInfo(changedFaces, changedFacesInfo);

    // Iterate until nothing changes
    iterate(maxIter);

    if ((maxIter > 0) && (iter_ >= maxIter))
    {
        FatalErrorIn
        (
            "FaceCellWave<Type>::FaceCellWave"
            "(const polyMesh&, const labelList&, const List<Type>,"
            " UList<Type>&, UList<Type>&, const label maxIter)"
        )
            << "Maximum number of iterations reached. Increase maxIter." << endl
            << "    maxIter:" << maxIter << endl
            << "    nChangedCells:" << nChangedCells_ << endl
            << "    nChangedFaces:" << nChangedFaces_ << endl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template <class Type>
Foam::label Foam::FaceCellWave<Type>::getUnsetCells() const
{
    return nUnvisitedCells_;
}


template <class Type>
Foam::label Foam::FaceCellWave<Type>::getUnsetFaces() const
{
    return nUnvisitedFaces_;
}



// Propagate cell to face
template <class Type>
Foam::label Foam::FaceCellWave<Type>::faceToCell()
{
    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();
    label nInternalFaces = mesh_.nInternalFaces();

    for
    (
        label changedFaceI = 0;
        changedFaceI < nChangedFaces_;
        changedFaceI++
    )
    {
        label faceI = changedFaces_[changedFaceI];
        if (!changedFace_[faceI])
        {
            FatalErrorIn("FaceCellWave<Type>::faceToCell()")
                << "Face " << faceI
                << " not marked as having been changed"
                << abort(FatalError);
        }


        const Type& neighbourWallInfo = allFaceInfo_[faceI];

        // Evaluate all connected cells

        // Owner
        label cellI = owner[faceI];
        Type& currentWallInfo = allCellInfo_[cellI];

        if (currentWallInfo != neighbourWallInfo)
        {
            updateCell
            (
                cellI,
                faceI,
                neighbourWallInfo,
                propagationTol_,
                currentWallInfo
            );
        }

        // Neighbour. Hack for check if face has neighbour.
        if (faceI < nInternalFaces)
        {
            cellI = neighbour[faceI];
            Type& currentWallInfo2 = allCellInfo_[cellI];

            if (currentWallInfo2 != neighbourWallInfo)
            {
                updateCell
                (
                    cellI,
                    faceI,
                    neighbourWallInfo,
                    propagationTol_,
                    currentWallInfo2
                );
            }
        }

        // Reset status of face
        changedFace_[faceI] = false;
    }

    // Handled all changed faces by now
    nChangedFaces_ = 0;

    if (debug)
    {
        Pout<< " Changed cells            : " << nChangedCells_ << endl;
    }

    // Sum nChangedCells over all procs
    label totNChanged = nChangedCells_;

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


// Propagate cell to face
template <class Type>
Foam::label Foam::FaceCellWave<Type>::cellToFace()
{
    const cellList& cells = mesh_.cells();

    for
    (
        label changedCellI = 0;
        changedCellI < nChangedCells_;
        changedCellI++
    )
    {
        label cellI = changedCells_[changedCellI];
        if (!changedCell_[cellI])
        {
            FatalErrorIn("FaceCellWave<Type>::cellToFace()") << "Cell " << cellI
                << " not marked as having been changed"
                << abort(FatalError);
        }

        const Type& neighbourWallInfo = allCellInfo_[cellI];

        // Evaluate all connected faces

        const labelList& faceLabels = cells[cellI];
        forAll(faceLabels, faceLabelI)
        {
            label faceI = faceLabels[faceLabelI];
            Type& currentWallInfo = allFaceInfo_[faceI];

            if (currentWallInfo != neighbourWallInfo)
            {
                updateFace
                (
                    faceI,
                    cellI,
                    neighbourWallInfo,
                    propagationTol_,
                    currentWallInfo
                );
            }
        }

        // Reset status of cell
        changedCell_[cellI] = false;
    }

    // Handled all changed cells by now
    nChangedCells_ = 0;

    if (hasCyclicPatches_)
    {
        // Transfer changed faces across cyclic halves
        handleCyclicPatches();
    }
    if (Pstream::parRun())
    {
        // Transfer changed faces from neighbouring processors.
        handleProcPatches();
    }

    if (debug)
    {
        Pout<< " Changed faces            : " << nChangedFaces_ << endl;
    }

    // Sum nChangedFaces over all procs
    label totNChanged = nChangedFaces_;

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


// Iterate
template <class Type>
Foam::label Foam::FaceCellWave<Type>::iterate(const label maxIter)
{
    if (hasCyclicPatches_)
    {
        // Transfer changed faces across cyclic halves
        handleCyclicPatches();
    }
    if (Pstream::parRun())
    {
        // Transfer changed faces from neighbouring processors.
        handleProcPatches();
    }

    while(iter_ < maxIter)
    {
        if (debug)
        {
            Pout<< " Iteration " << iter_ << endl;
        }

        nEvals_ = 0;

        label nCells = faceToCell();

        if (debug)
        {
            Pout<< " Total changed cells      : " << nCells << endl;
        }

        if (nCells == 0)
        {
            break;
        }

        label nFaces = cellToFace();

        if (debug)
        {
            Pout<< " Total changed faces      : " << nFaces << nl
                << " Total evaluations        : " << nEvals_ << nl
                << " Remaining unvisited cells: " << nUnvisitedCells_ << nl
                << " Remaining unvisited faces: " << nUnvisitedFaces_ << endl;
        }

        if (nFaces == 0)
        {
            break;
        }

        ++iter_;
    }

    return nUnvisitedCells_;
}

// ************************************************************************* //
