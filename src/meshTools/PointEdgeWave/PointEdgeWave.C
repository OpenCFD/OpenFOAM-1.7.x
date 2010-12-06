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

\*---------------------------------------------------------------------------*/

#include "PointEdgeWave.H"
#include "polyMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "OPstream.H"
#include "IPstream.H"
#include "PstreamCombineReduceOps.H"
#include "debug.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template <class Type>
Foam::scalar Foam::PointEdgeWave<Type>::propagationTol_ = 0.01;


// Offset labelList. Used for transferring from one cyclic half to the other.
template <class Type>
void Foam::PointEdgeWave<Type>::offset(const label val, labelList& elems)
{
    forAll(elems, i)
    {
        elems[i] += val;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Gets point-point correspondence. Is 
// - list of halfA points (in cyclic patch points)
// - list of halfB points (can overlap with A!)
// - for every patchPoint its corresponding point
template <class Type>
void Foam::PointEdgeWave<Type>::calcCyclicAddressing()
{
    label cycHalf = 0;

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];

        if (isA<cyclicPolyPatch>(patch))
        {
            label halfSize = patch.size()/2;

            SubList<face> halfAFaces
            (
                mesh_.faces(),
                halfSize,
                patch.start()
            );

            cycHalves_.set
            (
                cycHalf++,
                new primitivePatch(halfAFaces, mesh_.points())
            );

            SubList<face> halfBFaces
            (
                mesh_.faces(),
                halfSize,
                patch.start() + halfSize
            );

            cycHalves_.set
            (
                cycHalf++,
                new primitivePatch(halfBFaces, mesh_.points())
            );
        }
    }
}


// Handle leaving domain. Implementation referred to Type
template <class Type>
void Foam::PointEdgeWave<Type>::leaveDomain
(
    const polyPatch& meshPatch,
    const primitivePatch& patch,
    const labelList& patchPointLabels,
    List<Type>& pointInfo
) const
{
    const labelList& meshPoints = patch.meshPoints();
    
    forAll(patchPointLabels, i)
    {
        label patchPointI = patchPointLabels[i];

        const point& pt = patch.points()[meshPoints[patchPointI]];

        pointInfo[i].leaveDomain(meshPatch, patchPointI, pt);
    }
}


// Handle entering domain. Implementation referred to Type
template <class Type>
void Foam::PointEdgeWave<Type>::enterDomain
(
    const polyPatch& meshPatch,
    const primitivePatch& patch,
    const labelList& patchPointLabels,
    List<Type>& pointInfo
) const
{
    const labelList& meshPoints = patch.meshPoints();
    
    forAll(patchPointLabels, i)
    {
        label patchPointI = patchPointLabels[i];

        const point& pt = patch.points()[meshPoints[patchPointI]];

        pointInfo[i].enterDomain(meshPatch, patchPointI, pt);
    }
}


// Transform. Implementation referred to Type
template <class Type>
void Foam::PointEdgeWave<Type>::transform
(
    const tensorField& rotTensor,
    List<Type>& pointInfo
) const
{
    if (rotTensor.size() == 1)
    {
        const tensor& T = rotTensor[0];

        forAll(pointInfo, i)
        {
            pointInfo[i].transform(T);
        }
    }
    else
    {
        FatalErrorIn
        (
            "PointEdgeWave<Type>::transform(const tensorField&, List<Type>&)"
        )   << "Parallel cyclics not supported" << abort(FatalError);
    
        forAll(pointInfo, i)
        {
            pointInfo[i].transform(rotTensor[i]);
        }
    }
}


// Update info for pointI, at position pt, with information from
// neighbouring edge.
// Updates:
//      - changedPoint_, changedPoints_, nChangedPoints_,
//      - statistics: nEvals_, nUnvisitedPoints_
template <class Type>
bool Foam::PointEdgeWave<Type>::updatePoint
(
    const label pointI,
    const label neighbourEdgeI,
    const Type& neighbourInfo,
    const scalar tol,
    Type& pointInfo
)
{
    nEvals_++;

    bool wasValid = pointInfo.valid();

    bool propagate =
        pointInfo.updatePoint
        (
            mesh_,
            pointI,
            neighbourEdgeI,
            neighbourInfo,
            tol
        );

    if (propagate)
    {
        if (!changedPoint_[pointI])
        {
            changedPoint_[pointI] = true;
            changedPoints_[nChangedPoints_++] = pointI;
        }
    }

    if (!wasValid && pointInfo.valid())
    {
        --nUnvisitedPoints_;
    }

    return propagate;
}


// Update info for pointI, at position pt, with information from
// same point.
// Updates:
//      - changedPoint_, changedPoints_, nChangedPoints_,
//      - statistics: nEvals_, nUnvisitedPoints_
template <class Type>
bool Foam::PointEdgeWave<Type>::updatePoint
(
    const label pointI,
    const Type& neighbourInfo,
    const scalar tol,
    Type& pointInfo
)
{
    nEvals_++;

    bool wasValid = pointInfo.valid();

    bool propagate =
        pointInfo.updatePoint
        (
            mesh_,
            pointI,
            neighbourInfo,
            tol
        );

    if (propagate)
    {
        if (!changedPoint_[pointI])
        {
            changedPoint_[pointI] = true;
            changedPoints_[nChangedPoints_++] = pointI;
        }
    }

    if (!wasValid && pointInfo.valid())
    {
        --nUnvisitedPoints_;
    }

    return propagate;
}


// Update info for edgeI, at position pt, with information from
// neighbouring point.
// Updates:
//      - changedEdge_, changedEdges_, nChangedEdges_,
//      - statistics: nEvals_, nUnvisitedEdge_
template <class Type>
bool Foam::PointEdgeWave<Type>::updateEdge
(
    const label edgeI,
    const label neighbourPointI,
    const Type& neighbourInfo,
    const scalar tol,
    Type& edgeInfo
)
{
    nEvals_++;

    bool wasValid = edgeInfo.valid();

    bool propagate =
        edgeInfo.updateEdge
        (
            mesh_,
            edgeI,
            neighbourPointI,
            neighbourInfo,
            tol
        );

    if (propagate)
    {
        if (!changedEdge_[edgeI])
        {
            changedEdge_[edgeI] = true;
            changedEdges_[nChangedEdges_++] = edgeI;
        }
    }

    if (!wasValid && edgeInfo.valid())
    {
        --nUnvisitedEdges_;
    }

    return propagate;
}


// Check if patches of given type name are present
template <class Type>
template <class PatchType>
Foam::label Foam::PointEdgeWave<Type>::countPatchType() const
{
    label nPatches = 0;

    forAll(mesh_.boundaryMesh(), patchI)
    {
        if (isA<PatchType>(mesh_.boundaryMesh()[patchI]))
        {
            nPatches++;
        }
    }
    return nPatches;
}


// Collect changed patch points
template <class Type>
void Foam::PointEdgeWave<Type>::getChangedPatchPoints
(
    const primitivePatch& patch,

    DynamicList<Type>& patchInfo,
    DynamicList<label>& patchPoints,
    DynamicList<label>& owner,
    DynamicList<label>& ownerIndex
) const
{
    const labelList& meshPoints = patch.meshPoints();
    const faceList& localFaces = patch.localFaces();
    const labelListList& pointFaces = patch.pointFaces();

    forAll(meshPoints, patchPointI)
    {
        label meshPointI = meshPoints[patchPointI];

        if (changedPoint_[meshPointI])
        {
            patchInfo.append(allPointInfo_[meshPointI]);
            patchPoints.append(patchPointI);

            label patchFaceI = pointFaces[patchPointI][0];

            const face& f = localFaces[patchFaceI];

            label index = findIndex(f, patchPointI);

            owner.append(patchFaceI);
            ownerIndex.append(index);
        }
    }

    patchInfo.shrink();
    patchPoints.shrink();
    owner.shrink();
    ownerIndex.shrink();
}


// Update overall for changed patch points
template <class Type>
void Foam::PointEdgeWave<Type>::updateFromPatchInfo
(
    const polyPatch& meshPatch,
    const primitivePatch& patch,
    const labelList& owner,
    const labelList& ownerIndex,
    List<Type>& patchInfo
)
{
    const faceList& localFaces = patch.localFaces();
    const labelList& meshPoints = patch.meshPoints();

    // Get patch and mesh points.
    labelList changedPatchPoints(patchInfo.size());
    labelList changedMeshPoints(patchInfo.size());

    forAll(owner, i)
    {
        label faceI = owner[i];

        const face& f = localFaces[faceI];

        label index = (f.size() - ownerIndex[i]) % f.size();

        changedPatchPoints[i] = f[index];
        changedMeshPoints[i] = meshPoints[f[index]];
    }

    // Adapt for entering domain
    enterDomain(meshPatch, patch, changedPatchPoints, patchInfo);

    // Merge received info
    forAll(patchInfo, i)
    {
        updatePoint
        (
            changedMeshPoints[i],
            patchInfo[i],
            propagationTol_,
            allPointInfo_[changedMeshPoints[i]]
        );
    }
}



// Transfer all the information to/from neighbouring processors
template <class Type>
void Foam::PointEdgeWave<Type>::handleProcPatches()
{
    // 1. Send all point info on processor patches. Send as
    // face label + offset in face.

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(patch))
        {
            // Get all changed points in relative addressing

            DynamicList<Type> patchInfo(patch.nPoints());
            DynamicList<label> patchPoints(patch.nPoints());
            DynamicList<label> owner(patch.nPoints());
            DynamicList<label> ownerIndex(patch.nPoints());

            getChangedPatchPoints
            (
                patch,
                patchInfo,
                patchPoints,
                owner,
                ownerIndex
            );

            // Adapt for leaving domain
            leaveDomain(patch, patch, patchPoints, patchInfo);

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patch);

            if (debug)
            {
                Pout<< "Processor patch " << patchI << ' ' << patch.name()
                    << " communicating with " << procPatch.neighbProcNo()
                    << "  Sending:" << patchInfo.size() << endl;
            }

            {
                OPstream toNeighbour
                (   
                    Pstream::blocking,
                    procPatch.neighbProcNo()
                );

                toNeighbour << owner << ownerIndex << patchInfo;
            }
        }
    }


    //
    // 2. Receive all point info on processor patches.
    //

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];

        if (isA<processorPolyPatch>(patch))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patch);

            List<Type> patchInfo;
            labelList owner;
            labelList ownerIndex;
            {
                IPstream fromNeighbour
                (
                    Pstream::blocking,
                    procPatch.neighbProcNo()
                );

                fromNeighbour >> owner >> ownerIndex >> patchInfo;
            }    

            if (debug)
            {
                Pout<< "Processor patch " << patchI << ' ' << patch.name()
                    << " communicating with " << procPatch.neighbProcNo()
                    << "  Received:" << patchInfo.size() << endl;
            }

            // Apply transform to received data for non-parallel planes
            if (!procPatch.parallel())
            {
                transform(procPatch.forwardT(), patchInfo);
            }

            updateFromPatchInfo
            (
                patch,
                patch,
                owner,
                ownerIndex,
                patchInfo
            );
        }
    }



    //
    // 3. Handle all shared points
    //    (Note:irrespective if changed or not for now)
    //

    const globalMeshData& pd = mesh_.globalData();

    List<Type> sharedData(pd.nGlobalPoints());

    forAll(pd.sharedPointLabels(), i)
    {
        label meshPointI = pd.sharedPointLabels()[i];

        // Fill my entries in the shared points
        sharedData[pd.sharedPointAddr()[i]] = allPointInfo_[meshPointI];
    }

    // Combine on master. Reduce operator has to handle a list and call
    // Type.updatePoint for all elements
    combineReduce(sharedData, listUpdateOp<Type>());

    forAll(pd.sharedPointLabels(), i)
    {
        label meshPointI = pd.sharedPointLabels()[i];

        // Retrieve my entries from the shared points
        updatePoint
        (
            meshPointI,
            sharedData[pd.sharedPointAddr()[i]],
            propagationTol_,
            allPointInfo_[meshPointI]
        );
    }
}


template <class Type>
void Foam::PointEdgeWave<Type>::handleCyclicPatches()
{
    // 1. Send all point info on cyclic patches. Send as
    // face label + offset in face.

    label cycHalf = 0;

    forAll(mesh_.boundaryMesh(), patchI)
    {
        const polyPatch& patch = mesh_.boundaryMesh()[patchI];

        if (isA<cyclicPolyPatch>(patch))
        {
            const primitivePatch& halfA = cycHalves_[cycHalf++];
            const primitivePatch& halfB = cycHalves_[cycHalf++];

            // HalfA : get all changed points in relative addressing

            DynamicList<Type> halfAInfo(halfA.nPoints());
            DynamicList<label> halfAPoints(halfA.nPoints());
            DynamicList<label> halfAOwner(halfA.nPoints());
            DynamicList<label> halfAIndex(halfA.nPoints());

            getChangedPatchPoints
            (
                halfA,
                halfAInfo,
                halfAPoints,
                halfAOwner,
                halfAIndex
            );

            // HalfB : get all changed points in relative addressing

            DynamicList<Type> halfBInfo(halfB.nPoints());
            DynamicList<label> halfBPoints(halfB.nPoints());
            DynamicList<label> halfBOwner(halfB.nPoints());
            DynamicList<label> halfBIndex(halfB.nPoints());

            getChangedPatchPoints
            (
                halfB,
                halfBInfo,
                halfBPoints,
                halfBOwner,
                halfBIndex
            );


            // HalfA : adapt for leaving domain
            leaveDomain(patch, halfA, halfAPoints, halfAInfo);

            // HalfB : adapt for leaving domain
            leaveDomain(patch, halfB, halfBPoints, halfBInfo);


            // Apply rotation for non-parallel planes
            const cyclicPolyPatch& cycPatch =
                refCast<const cyclicPolyPatch>(patch);

            if (!cycPatch.parallel())
            {
                // received data from half1
                transform(cycPatch.reverseT(), halfAInfo);

                // received data from half2
                transform(cycPatch.forwardT(), halfBInfo);
            }

            if (debug)
            {
                Pout<< "Cyclic patch " << patchI << ' ' << patch.name() 
                    << "  Changed on first half : " << halfAInfo.size()
                    << "  Changed on second half : " << halfBInfo.size()
                    << endl;
            }

            // Half1: update with data from halfB
            updateFromPatchInfo
            (
                patch,
                halfA,
                halfBOwner,
                halfBIndex,
                halfBInfo
            );

            // Half2: update with data from halfA
            updateFromPatchInfo
            (
                patch,
                halfB,
                halfAOwner,
                halfAIndex,
                halfAInfo
            );

            if (debug)
            {
                //checkCyclic(patch);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Iterate, propagating changedPointsInfo across mesh, until no change (or 
// maxIter reached). Initial point values specified.
template <class Type>
Foam::PointEdgeWave<Type>::PointEdgeWave
(
    const polyMesh& mesh,
    const labelList& changedPoints,
    const List<Type>& changedPointsInfo,

    List<Type>& allPointInfo,
    List<Type>& allEdgeInfo,

    const label maxIter
)
:
    mesh_(mesh),
    allPointInfo_(allPointInfo),
    allEdgeInfo_(allEdgeInfo),
    changedPoint_(mesh_.nPoints(), false),
    changedPoints_(mesh_.nPoints()),
    nChangedPoints_(0),
    changedEdge_(mesh_.nEdges(), false),
    changedEdges_(mesh_.nEdges()),
    nChangedEdges_(0),
    nCyclicPatches_(countPatchType<cyclicPolyPatch>()),
    cycHalves_(2*nCyclicPatches_),
    nEvals_(0),
    nUnvisitedPoints_(mesh_.nPoints()),
    nUnvisitedEdges_(mesh_.nEdges())
{
    if (allPointInfo_.size() != mesh_.nPoints())
    {
        FatalErrorIn
        (
            "PointEdgeWave<Type>::PointEdgeWave"
            "(const polyMesh&, const labelList&, const List<Type>,"
            " List<Type>&, List<Type>&, const label maxIter)"
        )   << "size of pointInfo work array is not equal to the number"
            << " of points in the mesh" << endl
            << "    pointInfo   :" << allPointInfo_.size() << endl
            << "    mesh.nPoints:" << mesh_.nPoints()
            << exit(FatalError);
    }
    if (allEdgeInfo_.size() != mesh_.nEdges())
    {
        FatalErrorIn
        (
            "PointEdgeWave<Type>::PointEdgeWave"
            "(const polyMesh&, const labelList&, const List<Type>,"
            " List<Type>&, List<Type>&, const label maxIter)"
        )   << "size of edgeInfo work array is not equal to the number"
            << " of edges in the mesh" << endl
            << "    edgeInfo   :" << allEdgeInfo_.size() << endl
            << "    mesh.nEdges:" << mesh_.nEdges()
            << exit(FatalError);
    }


    // Calculate cyclic halves addressing.
    if (nCyclicPatches_ > 0)
    {
        calcCyclicAddressing();
    }


    // Set from initial changed points data
    setPointInfo(changedPoints, changedPointsInfo);

    if (debug)
    {
        Pout<< "Seed points               : " << nChangedPoints_ << endl;
    }

    // Iterate until nothing changes
    label iter = iterate(maxIter);

    if ((maxIter > 0) && (iter >= maxIter))
    {
        FatalErrorIn
        (
            "PointEdgeWave<Type>::PointEdgeWave"
            "(const polyMesh&, const labelList&, const List<Type>,"
            " List<Type>&, List<Type>&, const label maxIter)"
        )   << "Maximum number of iterations reached. Increase maxIter." << endl
            << "    maxIter:" << maxIter << endl
            << "    nChangedPoints:" << nChangedPoints_ << endl
            << "    nChangedEdges:" << nChangedEdges_ << endl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class Type>
Foam::PointEdgeWave<Type>::~PointEdgeWave()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template <class Type>
Foam::label Foam::PointEdgeWave<Type>::getUnsetPoints() const
{
    return nUnvisitedPoints_;
}


template <class Type>
Foam::label Foam::PointEdgeWave<Type>::getUnsetEdges() const
{
    return nUnvisitedEdges_;
}


// Copy point information into member data
template <class Type>
void Foam::PointEdgeWave<Type>::setPointInfo
(
    const labelList& changedPoints,
    const List<Type>& changedPointsInfo
)
{
    forAll(changedPoints, changedPointI)
    {
        label pointI = changedPoints[changedPointI];

        bool wasValid = allPointInfo_[pointI].valid();

        // Copy info for pointI
        allPointInfo_[pointI] = changedPointsInfo[changedPointI];

        // Maintain count of unset points
        if (!wasValid && allPointInfo_[pointI].valid())
        {
            --nUnvisitedPoints_;
        }

        // Mark pointI as changed, both on list and on point itself.

        if (!changedPoint_[pointI])
        {
            changedPoint_[pointI] = true;
            changedPoints_[nChangedPoints_++] = pointI;
        }
    }
}


// Propagate information from edge to point. Return number of points changed.
template <class Type>
Foam::label Foam::PointEdgeWave<Type>::edgeToPoint()
{
    for
    (
        label changedEdgeI = 0;
        changedEdgeI < nChangedEdges_;
        changedEdgeI++
    )
    {
        label edgeI = changedEdges_[changedEdgeI];

        if (!changedEdge_[edgeI])
        {
            FatalErrorIn("PointEdgeWave<Type>::edgeToPoint()")
                << "edge " << edgeI
                << " not marked as having been changed" << nl
                << "This might be caused by multiple occurences of the same"
                << " seed point." << abort(FatalError);
        }


        const Type& neighbourWallInfo = allEdgeInfo_[edgeI];

        // Evaluate all connected points (= edge endpoints)
        const edge& e = mesh_.edges()[edgeI];

        forAll(e, eI)
        {
            Type& currentWallInfo = allPointInfo_[e[eI]];

            if (currentWallInfo != neighbourWallInfo)
            {
                updatePoint
                (
                    e[eI],
                    edgeI,
                    neighbourWallInfo,
                    propagationTol_,
                    currentWallInfo
                );
            }
        }

        // Reset status of edge
        changedEdge_[edgeI] = false;
    }

    // Handled all changed edges by now
    nChangedEdges_ = 0;

    if (nCyclicPatches_ > 0)
    {
        // Transfer changed points across cyclic halves
        handleCyclicPatches();
    }
    if (Pstream::parRun())
    {
        // Transfer changed points from neighbouring processors.
        handleProcPatches();
    }

    if (debug)
    {
        Pout<< "Changed points            : " << nChangedPoints_ << endl;
    }

    // Sum nChangedPoints over all procs
    label totNChanged = nChangedPoints_;

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


// Propagate information from point to edge. Return number of edges changed.
template <class Type>
Foam::label Foam::PointEdgeWave<Type>::pointToEdge()
{
    const labelListList& pointEdges = mesh_.pointEdges();

    for
    (
        label changedPointI = 0;
        changedPointI < nChangedPoints_;
        changedPointI++
    )
    {
        label pointI = changedPoints_[changedPointI];

        if (!changedPoint_[pointI])
        {
            FatalErrorIn("PointEdgeWave<Type>::pointToEdge()")
                << "Point " << pointI
                << " not marked as having been changed" << nl
                << "This might be caused by multiple occurences of the same"
                << " seed point." << abort(FatalError);
        }

        const Type& neighbourWallInfo = allPointInfo_[pointI];

        // Evaluate all connected edges

        const labelList& edgeLabels = pointEdges[pointI];
        forAll(edgeLabels, edgeLabelI)
        {
            label edgeI = edgeLabels[edgeLabelI];

            Type& currentWallInfo = allEdgeInfo_[edgeI];

            if (currentWallInfo != neighbourWallInfo)
            {
                updateEdge
                (
                    edgeI,
                    pointI,
                    neighbourWallInfo,
                    propagationTol_,
                    currentWallInfo
                );
            }
        }

        // Reset status of point
        changedPoint_[pointI] = false;
    }

    // Handled all changed points by now
    nChangedPoints_ = 0;

    if (debug)
    {
        Pout<< "Changed edges             : " << nChangedEdges_ << endl;
    }

    // Sum nChangedPoints over all procs
    label totNChanged = nChangedEdges_;

    reduce(totNChanged, sumOp<label>());

    return totNChanged;
}


// Iterate
template <class Type>
Foam::label Foam::PointEdgeWave<Type>::iterate(const label maxIter)
{
    if (nCyclicPatches_ > 0)
    {
        // Transfer changed points across cyclic halves
        handleCyclicPatches();
    }
    if (Pstream::parRun())
    {
        // Transfer changed points from neighbouring processors.
        handleProcPatches();
    }

    nEvals_ = 0;

    label iter = 0;

    while(iter < maxIter)
    {
        if (debug)
        {
            Pout<< "Iteration " << iter << endl;
        }

        label nEdges = pointToEdge();

        if (debug)
        {
            Pout<< "Total changed edges       : " << nEdges << endl;
        }

        if (nEdges == 0)
        {
            break;
        }

        label nPoints = edgeToPoint();

        if (debug)
        {
            Pout<< "Total changed points      : " << nPoints << endl;

            Pout<< "Total evaluations         : " << nEvals_ << endl;

            Pout<< "Remaining unvisited points: " << nUnvisitedPoints_ << endl;

            Pout<< "Remaining unvisited edges : " << nUnvisitedEdges_ << endl;

            Pout<< endl;
        }

        if (nPoints == 0)
        {
            break;
        }

        iter++;
    }

    return iter;
}

// ************************************************************************* //
