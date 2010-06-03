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

#include "globalMeshData.H"
#include "Time.H"
#include "Pstream.H"
#include "PstreamCombineReduceOps.H"
#include "processorPolyPatch.H"
#include "demandDrivenData.H"
#include "globalPoints.H"
//#include "geomGlobalPoints.H"
#include "labelIOList.H"
#include "PackedList.H"
#include "mergePoints.H"
#include "matchPoints.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::globalMeshData, 0);

// Geometric matching tolerance. Factor of mesh bounding box.
const Foam::scalar Foam::globalMeshData::matchTol_ = 1E-8;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Collect processor patch addressing.
void Foam::globalMeshData::initProcAddr()
{
    processorPatchIndices_.setSize(mesh_.boundaryMesh().size());
    processorPatchIndices_ = -1;

    processorPatchNeighbours_.setSize(mesh_.boundaryMesh().size());
    processorPatchNeighbours_ = -1;

    // Construct processor patch indexing. processorPatchNeighbours_ only
    // set if running in parallel!
    processorPatches_.setSize(mesh_.boundaryMesh().size());

    label nNeighbours = 0;

    forAll (mesh_.boundaryMesh(), patchi)
    {
        if (isA<processorPolyPatch>(mesh_.boundaryMesh()[patchi]))
        {
            processorPatches_[nNeighbours] = patchi;
            processorPatchIndices_[patchi] = nNeighbours++;
        }
    }
    processorPatches_.setSize(nNeighbours);


    if (Pstream::parRun())
    {
        // Send indices of my processor patches to my neighbours
        forAll (processorPatches_, i)
        {
            label patchi = processorPatches_[i];

            OPstream toNeighbour
            (
                Pstream::blocking,
                refCast<const processorPolyPatch>
                (
                    mesh_.boundaryMesh()[patchi]
                ).neighbProcNo()
            );

            toNeighbour << processorPatchIndices_[patchi];
        }

        forAll(processorPatches_, i)
        {
            label patchi = processorPatches_[i];
            
            IPstream fromNeighbour
            (
                Pstream::blocking,
                refCast<const processorPolyPatch>
                (
                    mesh_.boundaryMesh()[patchi]
                ).neighbProcNo()
            );
            
            fromNeighbour >> processorPatchNeighbours_[patchi];
        }
    }
}


// Given information about locally used edges allocate global shared edges.
void Foam::globalMeshData::countSharedEdges
(
    const HashTable<labelList, edge, Hash<edge> >& procSharedEdges,
    HashTable<label, edge, Hash<edge> >& globalShared,
    label& sharedEdgeI
)
{
    // Count occurrences of procSharedEdges in global shared edges table.
    for
    (
        HashTable<labelList, edge, Hash<edge> >::const_iterator iter =
            procSharedEdges.begin();
        iter != procSharedEdges.end();
        ++iter
    )
    {
        const edge& e = iter.key();

        HashTable<label, edge, Hash<edge> >::iterator globalFnd =
            globalShared.find(e);

        if (globalFnd == globalShared.end())
        {
            // First time occurrence of this edge. Check how many we are adding.
            if (iter().size() == 1)
            {
                // Only one edge. Mark with special value.
                globalShared.insert(e, -1);
            }
            else
            {
                // Edge used more than once (even by local shared edges alone)
                // so allocate proper shared edge label.
                globalShared.insert(e, sharedEdgeI++);
            }
        }
        else
        {
            if (globalFnd() == -1)
            {
                // Second time occurence of this edge. Assign proper
                // edge label.
                globalFnd() = sharedEdgeI++;
            }
        }
    }
}


// Shared edges are shared between multiple processors. By their nature both
// of their endpoints are shared points. (but not all edges using two shared
// points are shared edges! There might e.g. be an edge between two unrelated
// clusters of shared points)
void Foam::globalMeshData::calcSharedEdges() const
{
    if (nGlobalEdges_ != -1 || sharedEdgeLabelsPtr_ || sharedEdgeAddrPtr_)
    {
        FatalErrorIn("globalMeshData::calcSharedEdges()")
            << "Shared edge addressing already done" << abort(FatalError);
    }


    const labelList& sharedPtAddr = sharedPointAddr();
    const labelList& sharedPtLabels = sharedPointLabels();

    // Since don't want to construct pointEdges for whole mesh create
    // Map for all shared points.
    Map<label> meshToShared(2*sharedPtLabels.size());
    forAll(sharedPtLabels, i)
    {
        meshToShared.insert(sharedPtLabels[i], i);
    }

    // Find edges using shared points. Store correspondence to local edge
    // numbering. Note that multiple local edges can have the same shared
    // points! (for cyclics or separated processor patches)
    HashTable<labelList, edge, Hash<edge> > localShared
    (
        2*sharedPtAddr.size()
    );

    const edgeList& edges = mesh_.edges();

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        Map<label>::const_iterator e0Fnd = meshToShared.find(e[0]);

        if (e0Fnd != meshToShared.end())
        {
            Map<label>::const_iterator e1Fnd = meshToShared.find(e[1]);

            if (e1Fnd != meshToShared.end())
            {
                // Found edge which uses shared points. Probably shared.

                // Construct the edge in shared points (or rather global indices
                // of the shared points)
                edge sharedEdge
                (
                    sharedPtAddr[e0Fnd()],
                    sharedPtAddr[e1Fnd()]
                );

                HashTable<labelList, edge, Hash<edge> >::iterator iter =
                    localShared.find(sharedEdge);

                if (iter == localShared.end())
                {
                    // First occurrence of this point combination. Store.
                    localShared.insert(sharedEdge, labelList(1, edgeI));
                }
                else
                {
                    // Add this edge to list of edge labels.
                    labelList& edgeLabels = iter();

                    label sz = edgeLabels.size();
                    edgeLabels.setSize(sz+1);
                    edgeLabels[sz] = edgeI;
                }
            }
        }
    }


    // Now we have a table on every processors which gives its edges which use
    // shared points. Send this all to the master and have it allocate
    // global edge numbers for it. But only allocate a global edge number for
    // edge if it is used more than once!
    // Note that we are now sending the whole localShared to the master whereas
    // we only need the local count (i.e. the number of times a global edge is
    // used). But then this only gets done once so not too bothered about the
    // extra global communication.

    HashTable<label, edge, Hash<edge> > globalShared(nGlobalPoints());

    if (Pstream::master())
    {
        label sharedEdgeI = 0;

        // Merge my shared edges into the global list
        if (debug)
        {
            Pout<< "globalMeshData::calcSharedEdges : Merging in from proc0 : "
                << localShared.size() << endl;
        }
        countSharedEdges(localShared, globalShared, sharedEdgeI);

        // Receive data from slaves and insert
        if (Pstream::parRun())
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                // Receive the edges using shared points from the slave.
                IPstream fromSlave(Pstream::blocking, slave);
                HashTable<labelList, edge, Hash<edge> > procSharedEdges
                (
                    fromSlave
                );

                if (debug)
                {
                    Pout<< "globalMeshData::calcSharedEdges : "
                        << "Merging in from proc"
                        << Foam::name(slave) << " : " << procSharedEdges.size()
                        << endl;
                }
                countSharedEdges(procSharedEdges, globalShared, sharedEdgeI);
            }
        }

        // Now our globalShared should have some edges with -1 as edge label
        // These were only used once so are not proper shared edges.
        // Remove them.
        {
            HashTable<label, edge, Hash<edge> > oldSharedEdges(globalShared);

            globalShared.clear();

            for
            (
                HashTable<label, edge, Hash<edge> >::const_iterator iter =
                    oldSharedEdges.begin();
                iter != oldSharedEdges.end();
                ++iter
            )
            {
                if (iter() != -1)
                {
                    globalShared.insert(iter.key(), iter());
                }
            }
            if (debug)
            {
                Pout<< "globalMeshData::calcSharedEdges : Filtered "
                    << oldSharedEdges.size()
                    << " down to " << globalShared.size() << endl;
            }
        }


        // Send back to slaves.
        if (Pstream::parRun())
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                // Receive the edges using shared points from the slave.
                OPstream toSlave(Pstream::blocking, slave);
                toSlave << globalShared;
            }
        }
    }
    else
    {
        // Send local edges to master
        {
            OPstream toMaster(Pstream::blocking, Pstream::masterNo());

            toMaster << localShared;
        }
        // Receive merged edges from master.
        {
            IPstream fromMaster(Pstream::blocking, Pstream::masterNo());

            fromMaster >> globalShared;
        }
    }

    // Now use the global shared edges list (globalShared) to classify my local
    // ones (localShared)

    nGlobalEdges_ = globalShared.size();

    DynamicList<label> dynSharedEdgeLabels(globalShared.size());
    DynamicList<label> dynSharedEdgeAddr(globalShared.size());

    for
    (
        HashTable<labelList, edge, Hash<edge> >::const_iterator iter =
            localShared.begin();
        iter != localShared.end();
        ++iter
    )
    {
        const edge& e = iter.key();

        HashTable<label, edge, Hash<edge> >::const_iterator edgeFnd =
            globalShared.find(e);

        if (edgeFnd != globalShared.end())
        {
            // My local edge is indeed a shared one. Go through all local edge
            // labels with this point combination.
            const labelList& edgeLabels = iter();

            forAll(edgeLabels, i)
            {
                // Store label of local mesh edge
                dynSharedEdgeLabels.append(edgeLabels[i]);

                // Store label of shared edge
                dynSharedEdgeAddr.append(edgeFnd());
            }
        }
    }

    sharedEdgeLabelsPtr_ = new labelList();
    labelList& sharedEdgeLabels = *sharedEdgeLabelsPtr_;
    sharedEdgeLabels.transfer(dynSharedEdgeLabels);

    sharedEdgeAddrPtr_ = new labelList();
    labelList& sharedEdgeAddr = *sharedEdgeAddrPtr_;
    sharedEdgeAddr.transfer(dynSharedEdgeAddr);

    if (debug)
    {
        Pout<< "globalMeshData : nGlobalEdges_:" << nGlobalEdges_ << nl
            << "globalMeshData : sharedEdgeLabels:" << sharedEdgeLabels.size()
            << nl
            << "globalMeshData : sharedEdgeAddr:" << sharedEdgeAddr.size()
            << endl;
    }
}


// Helper function to count coincident faces. This part used to be
// in updateMesh but I've had to move it to a separate function
// because of aliasing optimisation errors in icc9.1 on the
// Itanium.
Foam::label Foam::globalMeshData::countCoincidentFaces
(
    const scalar tolDim,
    const vectorField& separationDist
)
{
    label nCoincident = 0;

    forAll(separationDist, faceI)
    {
        if (mag(separationDist[faceI]) < tolDim)
        {
            // Faces coincide
            nCoincident++;
        }
    }
    return nCoincident;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from polyMesh
Foam::globalMeshData::globalMeshData(const polyMesh& mesh)
:
    processorTopology(mesh.boundaryMesh()),
    mesh_(mesh),
    bb_(vector::zero, vector::zero),
    nTotalPoints_(-1),
    nTotalFaces_(-1),
    nTotalCells_(-1),
    processorPatches_(0),
    processorPatchIndices_(0),
    processorPatchNeighbours_(0),
    nGlobalPoints_(-1),
    sharedPointLabels_(0),
    sharedPointAddr_(0),
    sharedPointGlobalLabelsPtr_(NULL),
    nGlobalEdges_(-1),
    sharedEdgeLabelsPtr_(NULL),
    sharedEdgeAddrPtr_(NULL)
{
    updateMesh();
}


// Read constructor given IOobject and a polyMesh reference
Foam::globalMeshData::globalMeshData(const IOobject& io, const polyMesh& mesh)
:
    processorTopology(mesh.boundaryMesh()),
    mesh_(mesh),
    bb_(mesh.points()),
    nTotalPoints_(-1),
    nTotalFaces_(-1),
    nTotalCells_(-1),
    processorPatches_(0),
    processorPatchIndices_(0),
    processorPatchNeighbours_(0),
    nGlobalPoints_(-1),
    sharedPointLabels_(0),
    sharedPointAddr_(0),
    sharedPointGlobalLabelsPtr_(NULL),
    nGlobalEdges_(-1),
    sharedEdgeLabelsPtr_(NULL),
    sharedEdgeAddrPtr_(NULL)
{
    initProcAddr();

    IOdictionary dict(io);

    dict.lookup("nTotalPoints") >> nTotalPoints_;
    dict.lookup("nTotalFaces") >> nTotalFaces_;
    dict.lookup("nTotalCells") >> nTotalCells_;
    dict.lookup("nGlobalPoints") >> nGlobalPoints_;
    dict.lookup("sharedPointLabels") >> sharedPointLabels_;
    dict.lookup("sharedPointAddr") >> sharedPointAddr_;
    labelList sharedPointGlobalLabels(dict.lookup("sharedPointGlobalLabels"));

    sharedPointGlobalLabelsPtr_ = new labelList(sharedPointGlobalLabels);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::globalMeshData::~globalMeshData()
{
    clearOut();
}


void Foam::globalMeshData::clearOut()
{
    deleteDemandDrivenData(sharedPointGlobalLabelsPtr_);
    // Edge
    nGlobalPoints_ = -1;
    deleteDemandDrivenData(sharedEdgeLabelsPtr_);
    deleteDemandDrivenData(sharedEdgeAddrPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return shared point global labels.
const Foam::labelList& Foam::globalMeshData::sharedPointGlobalLabels() const
{
    if (!sharedPointGlobalLabelsPtr_)
    {
        sharedPointGlobalLabelsPtr_ = new labelList(sharedPointLabels_.size());
        labelList& sharedPointGlobalLabels = *sharedPointGlobalLabelsPtr_;

        IOobject addrHeader
        (
            "pointProcAddressing",
            mesh_.facesInstance()/mesh_.meshSubDir,
            mesh_,
            IOobject::MUST_READ
        );

        if (addrHeader.headerOk())
        {
            // There is a pointProcAddressing file so use it to get labels
            // on the original mesh
            Pout<< "globalMeshData::sharedPointGlobalLabels : "
                << "Reading pointProcAddressing" << endl;

            labelIOList pointProcAddressing(addrHeader);

            forAll(sharedPointLabels_, i)
            {
                // Get my mesh point
                label pointI = sharedPointLabels_[i];

                // Map to mesh point of original mesh
                sharedPointGlobalLabels[i] = pointProcAddressing[pointI];
            }
        }
        else
        {
            Pout<< "globalMeshData::sharedPointGlobalLabels :"
                << " Setting pointProcAddressing to -1" << endl;

            sharedPointGlobalLabels = -1;
        }
    }
    return *sharedPointGlobalLabelsPtr_;
}


// Collect coordinates of shared points. (does parallel communication!)
Foam::pointField Foam::globalMeshData::sharedPoints() const
{
    // Get all processors to send their shared points to master.
    // (not very efficient)

    pointField sharedPoints(nGlobalPoints_);

    if (Pstream::master())
    {
        // Master:
        // insert my own data first
        forAll(sharedPointLabels_, i)
        {
            label sharedPointI = sharedPointAddr_[i];

            sharedPoints[sharedPointI] = mesh_.points()[sharedPointLabels_[i]];
        }

        // Receive data from slaves and insert
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            IPstream fromSlave(Pstream::blocking, slave);

            labelList nbrSharedPointAddr;
            pointField nbrSharedPoints;
            fromSlave >> nbrSharedPointAddr >> nbrSharedPoints;

            forAll(nbrSharedPointAddr, i)
            {
                label sharedPointI = nbrSharedPointAddr[i];

                sharedPoints[sharedPointI] = nbrSharedPoints[i];
            }
        }

        // Send back
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            OPstream toSlave
            (
                Pstream::blocking,
                slave,
                sharedPoints.size()*sizeof(vector::zero)
            );
            toSlave << sharedPoints;
        }
    }
    else
    {
        // Slave:
        // send points
        {
            OPstream toMaster(Pstream::blocking, Pstream::masterNo());

            toMaster
                << sharedPointAddr_ 
                << UIndirectList<point>(mesh_.points(), sharedPointLabels_)();
        }

        // Receive sharedPoints
        {
            IPstream fromMaster(Pstream::blocking, Pstream::masterNo());
            fromMaster >> sharedPoints;
        }
    }

    return sharedPoints;
}


// Collect coordinates of shared points. (does parallel communication!)
Foam::pointField Foam::globalMeshData::geometricSharedPoints() const
{
    // Get coords of my shared points
    pointField sharedPoints(sharedPointLabels_.size());

    forAll(sharedPointLabels_, i)
    {
        label meshPointI = sharedPointLabels_[i];

        sharedPoints[i] = mesh_.points()[meshPointI];
    }

    // Append from all processors
    combineReduce(sharedPoints, plusEqOp<pointField>());

    // Merge tolerance
    scalar tolDim = matchTol_ * bb_.mag();

    // And see how many are unique
    labelList pMap;
    pointField mergedPoints;

    mergePoints
    (
        sharedPoints,   // coordinates to merge
        tolDim,         // tolerance
        false,          // verbosity
        pMap,
        mergedPoints
    );

    return mergedPoints;
}


Foam::label Foam::globalMeshData::nGlobalEdges() const
{
    if (nGlobalEdges_ == -1)
    {
        calcSharedEdges();
    }
    return nGlobalEdges_;
}


const Foam::labelList& Foam::globalMeshData::sharedEdgeLabels() const
{
    if (!sharedEdgeLabelsPtr_)
    {
        calcSharedEdges();
    }
    return *sharedEdgeLabelsPtr_;
}


const Foam::labelList& Foam::globalMeshData::sharedEdgeAddr() const
{
    if (!sharedEdgeAddrPtr_)
    {
        calcSharedEdges();
    }
    return *sharedEdgeAddrPtr_;
}


void Foam::globalMeshData::movePoints(const pointField& newPoints)
{
    // Topology does not change and we don't store any geometry so nothing
    // needs to be done.
}


// Update all data after morph
void Foam::globalMeshData::updateMesh()
{
    // Clear out old data
    clearOut();

    // Do processor patch addressing
    initProcAddr();

    // Note: boundBox does reduce
    bb_ = boundBox(mesh_.points());

    scalar tolDim = matchTol_ * bb_.mag();

    if (debug)
    {
        Pout<< "globalMeshData : bb_:" << bb_
            << " merge dist:" << tolDim << endl;
    }


    // Option 1. Topological
    {
        // Calculate all shared points. This does all the hard work.
        globalPoints parallelPoints(mesh_);

        // Copy data out.
        nGlobalPoints_ = parallelPoints.nGlobalPoints();
        sharedPointLabels_ = parallelPoints.sharedPointLabels();
        sharedPointAddr_ = parallelPoints.sharedPointAddr();
    }
    //// Option 2. Geometric
    //{
    //    // Calculate all shared points. This does all the hard work.
    //    geomGlobalPoints parallelPoints(mesh_, tolDim);
    //
    //    // Copy data out.
    //    nGlobalPoints_ = parallelPoints.nGlobalPoints();
    //    sharedPointLabels_ = parallelPoints.sharedPointLabels();
    //    sharedPointAddr_ = parallelPoints.sharedPointAddr();
    //
    //    nGlobalEdges_ = parallelPoints.nGlobalEdges();
    //    sharedEdgeLabels_ = parallelPoints.sharedEdgeLabels();
    //    sharedEdgeAddr_ = parallelPoints.sharedEdgeAddr();
    //}

    // Total number of faces. Start off from all faces. Remove coincident
    // processor faces (on highest numbered processor) before summing.
    nTotalFaces_ = mesh_.nFaces();

    // Do not count processor-patch faces that are coincident.
    forAll(processorPatches_, i)
    {
        label patchI = processorPatches_[i];

        const processorPolyPatch& procPatch =
            refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

        if (Pstream::myProcNo() > procPatch.neighbProcNo())
        {
            // Uncount my faces. Handle cyclics separately.

            if (procPatch.separated())
            {
                const vectorField& separationDist = procPatch.separation();

                nTotalFaces_ -= countCoincidentFaces(tolDim, separationDist);
            }
            else
            {
                // Normal, unseparated processor patch. Remove duplicates.
                nTotalFaces_ -= procPatch.size();
            }
        }
    }
    reduce(nTotalFaces_, sumOp<label>());

    if (debug)
    {
        Pout<< "globalMeshData : nTotalFaces_:" << nTotalFaces_ << endl;
    }


    nTotalCells_ = mesh_.nCells();
    reduce(nTotalCells_, sumOp<label>());

    if (debug)
    {
        Pout<< "globalMeshData : nTotalCells_:" << nTotalCells_ << endl;
    }

    nTotalPoints_ = mesh_.nPoints();

    // Correct points for duplicate ones. We have
    // - points shared between 2 processor patches only. Count only on
    //   lower numbered processor. Make sure to count only once since points
    //   can be on multiple patches on the same processor.
    // - globally shared points.

    if (Pstream::parRun())
    {
        const label UNSET = 0;      // not set
        const label SHARED = 1;     // globally shared
        const label VISITED = 2;    // corrected for

        // Mark globally shared points
        PackedList<2> pointStatus(mesh_.nPoints(), UNSET);

        forAll(sharedPointLabels_, i)
        {
            label meshPointI = sharedPointLabels_[i];

            pointStatus.set(meshPointI, SHARED);
        }

        // Send patch local points
        forAll(processorPatches_, i)
        {
            label patchI = processorPatches_[i];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

            OPstream toNeighbour(Pstream::blocking, procPatch.neighbProcNo());

            toNeighbour << procPatch.localPoints();
        }

        // Receive patch local points and uncount if coincident (and not shared)
        forAll(processorPatches_, i)
        {
            label patchI = processorPatches_[i];

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(mesh_.boundaryMesh()[patchI]);

            IPstream fromNeighbour(Pstream::blocking, procPatch.neighbProcNo());

            pointField nbrPoints(fromNeighbour);

            if (Pstream::myProcNo() > procPatch.neighbProcNo())
            {
                labelList pMap;
                matchPoints
                (
                    procPatch.localPoints(),
                    nbrPoints,
                    scalarField(procPatch.nPoints(), tolDim),   // tolerance
                    false,      // verbosity
                    pMap        // map from my points to nbrPoints
                );

                forAll(pMap, patchPointI)
                {
                    label meshPointI = procPatch.meshPoints()[patchPointI];

                    label stat = pointStatus.get(meshPointI);

                    if (stat == UNSET)
                    {
                        // Mark point as visited so if point is on multiple proc
                        // patches it only gets uncounted once.
                        pointStatus.set(meshPointI, VISITED);

                        if (pMap[patchPointI] != -1)
                        {
                            // Points share same coordinate so uncount.
                            nTotalPoints_--;
                        }
                    }
                }
            }
        }
        // Sum all points
        reduce(nTotalPoints_, sumOp<label>());
    }

    // nTotalPoints has not been corrected yet for shared points. For these
    // just collect all their coordinates and count unique ones.

    label mySharedPoints = sharedPointLabels_.size();
    reduce(mySharedPoints, sumOp<label>());

    // Collect and merge shared points (does parallel communication)
    pointField geomSharedPoints(geometricSharedPoints());
    label nGeomSharedPoints = geomSharedPoints.size();

    // Shared points merged down to mergedPoints size.
    nTotalPoints_ -= mySharedPoints - nGeomSharedPoints;

    if (debug)
    {
        Pout<< "globalMeshData : nTotalPoints_:" << nTotalPoints_ << endl;
    }

    //
    // Now we have all info we wanted.
    // Do some checking (if debug is set)
    //

    if (debug)
    {
        if (Pstream::master())
        {
            // We have the geometricSharedPoints already so write them.
            // Ideally would like to write the networks of connected points as
            // well but this is harder. (Todo)
            Pout<< "globalMeshData : writing geometrically separated shared"
                << " points to geomSharedPoints.obj" << endl;

            OFstream str("geomSharedPoints.obj");

            forAll(geomSharedPoints, i)
            {
                const point& pt = geomSharedPoints[i];

                str << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z()
                    << nl;
            }
        }
    }
}


// Write data
bool Foam::globalMeshData::write() const
{
    IOdictionary dict
    (
        IOobject
        (
            "parallelData",
            mesh_.facesInstance(),
            mesh_.meshSubDir,
            mesh_
        )
    );

    dict.add("nTotalPoints", nTotalPoints());
    dict.add("nTotalFaces", nTotalFaces());
    dict.add("nTotalCells", nTotalCells());

    dict.add("nGlobalPoints", nGlobalPoints());
    dict.add("sharedPointLabels", sharedPointLabels());
    dict.add("sharedPointAddr", sharedPointAddr());
    dict.add("sharedPointGlobalLabels", sharedPointGlobalLabels());

    return dict.writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED
    );
}


// * * * * * * * * * * * * * * * Ostream Operators * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const globalMeshData& p)
{
    os  << "nTotalPoints " << p.nTotalPoints() << token::END_STATEMENT << nl
        << "nTotalFaces " << p.nTotalFaces() << token::END_STATEMENT << nl
        << "nTotalCells " << p.nTotalCells() << token::END_STATEMENT << nl
        << "nGlobalPoints " << p.nGlobalPoints() << token::END_STATEMENT << nl
        << "sharedPointLabels " << p.sharedPointLabels()
        << token::END_STATEMENT << nl
        << "sharedPointAddr " << p.sharedPointAddr()
        << token::END_STATEMENT << endl;

    return os;
}


// ************************************************************************* //
