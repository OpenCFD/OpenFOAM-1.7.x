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

Description
    All to do with snapping to the surface

\*----------------------------------------------------------------------------*/

#include "autoSnapDriver.H"
#include "Time.H"
#include "motionSmoother.H"
#include "polyTopoChange.H"
#include "OFstream.H"
#include "syncTools.H"
#include "fvMesh.H"
#include "Time.H"
#include "OFstream.H"
#include "mapPolyMesh.H"
#include "motionSmoother.H"
#include "pointEdgePoint.H"
#include "PointEdgeWave.H"
#include "mergePoints.H"
#include "snapParameters.H"
#include "refinementSurfaces.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(autoSnapDriver, 0);

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Get faces to repatch. Returns map from face to patch.
Foam::Map<Foam::label> Foam::autoSnapDriver::getZoneBafflePatches
(
    const bool allowBoundary
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    Map<label> bafflePatch(mesh.nFaces()/1000);

    const wordList& faceZoneNames = surfaces.faceZoneNames();
    const faceZoneMesh& fZones = mesh.faceZones();

    forAll(faceZoneNames, surfI)
    {
        if (faceZoneNames[surfI].size())
        {
            // Get zone
            label zoneI = fZones.findZoneID(faceZoneNames[surfI]);

            const faceZone& fZone = fZones[zoneI];

            //// Get patch allocated for zone
            //label patchI = surfaceToCyclicPatch_[surfI];
            // Get patch of (first region) of surface
            label patchI = globalToPatch_[surfaces.globalRegion(surfI, 0)];

            Info<< "For surface "
                << surfaces.names()[surfI]
                << " found faceZone " << fZone.name()
                << " and patch " << mesh.boundaryMesh()[patchI].name()
                << endl;


            forAll(fZone, i)
            {
                label faceI = fZone[i];

                if (allowBoundary || mesh.isInternalFace(faceI))
                {
                    if (!bafflePatch.insert(faceI, patchI))
                    {
                        label oldPatchI = bafflePatch[faceI];

                        if (oldPatchI != patchI)
                        {
                            FatalErrorIn("getZoneBafflePatches(const bool)")
                                << "Face " << faceI
                                << " fc:" << mesh.faceCentres()[faceI]
                                << " is in faceZone "
                                << mesh.boundaryMesh()[oldPatchI].name()
                                << " and in faceZone "
                                << mesh.boundaryMesh()[patchI].name()
                                << abort(FatalError);
                        }
                    }
                }
            }
        }
    }
    return bafflePatch;
}


// Calculate geometrically collocated points, Requires PackedList to be
// sized and initalised!
Foam::label Foam::autoSnapDriver::getCollocatedPoints
(
    const scalar tol,
    const pointField& points,
    PackedBoolList& isCollocatedPoint
)
{
    labelList pointMap;
    pointField newPoints;
    bool hasMerged = mergePoints
    (
        points,                         // points
        tol,                            // mergeTol
        false,                          // verbose
        pointMap,
        newPoints
    );

    if (!returnReduce(hasMerged, orOp<bool>()))
    {
        return 0;
    }

    // Determine which newPoints are referenced more than once
    label nCollocated = 0;

    // Per old point the newPoint. Or -1 (not set yet) or -2 (already seen
    // twice)
    labelList firstOldPoint(newPoints.size(), -1);
    forAll(pointMap, oldPointI)
    {
        label newPointI = pointMap[oldPointI];

        if (firstOldPoint[newPointI] == -1)
        {
            // First use of oldPointI. Store.
            firstOldPoint[newPointI] = oldPointI;
        }
        else if (firstOldPoint[newPointI] == -2)
        {
            // Third or more reference of oldPointI -> non-manifold
            isCollocatedPoint.set(oldPointI, 1u);
            nCollocated++;
        }
        else
        {
            // Second reference of oldPointI -> non-manifold
            isCollocatedPoint.set(firstOldPoint[newPointI], 1u);
            nCollocated++;

            isCollocatedPoint.set(oldPointI, 1u);
            nCollocated++;

            // Mark with special value to save checking next time round
            firstOldPoint[newPointI] = -2;
        }
    }
    return returnReduce(nCollocated, sumOp<label>());
}


// Calculate displacement as average of patch points.
Foam::pointField Foam::autoSnapDriver::smoothPatchDisplacement
(
    const motionSmoother& meshMover,
    const List<labelPair>& baffles
) const
{
    const indirectPrimitivePatch& pp = meshMover.patch();

    // Calculate geometrically non-manifold points on the patch to be moved.
    PackedBoolList nonManifoldPoint(pp.nPoints());
    label nNonManifoldPoints = getCollocatedPoints
    (
        SMALL,
        pp.localPoints(),
        nonManifoldPoint
    );
    Info<< "Found " << nNonManifoldPoints << " non-mainfold point(s)."
        << endl;


    // Average points
    // ~~~~~~~~~~~~~~

    // We determine three points:
    // - average of (centres of) connected patch faces
    // - average of (centres of) connected internal mesh faces
    // - as fallback: centre of any connected cell
    // so we can do something moderately sensible for non/manifold points.

    // Note: the averages are calculated properly parallel. This is
    // necessary to get the points shared by processors correct.


    const labelListList& pointFaces = pp.pointFaces();
    const labelList& meshPoints = pp.meshPoints();
    const pointField& points = pp.points();
    const polyMesh& mesh = meshMover.mesh();

    // Get labels of faces to count (master of coupled faces and baffle pairs)
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));

    {
        forAll(baffles, i)
        {
            label f0 = baffles[i].first();
            label f1 = baffles[i].second();

            if (isMasterFace.get(f0) == 1)
            {
                // Make f1 a slave
                isMasterFace.set(f1, 0);
            }
            else if (isMasterFace.get(f1) == 1)
            {
                isMasterFace.set(f0, 0);
            }
            else
            {
                FatalErrorIn("autoSnapDriver::smoothPatchDisplacement(..)")
                    << "Both sides of baffle consisting of faces " << f0
                    << " and " << f1 << " are already slave faces."
                    << abort(FatalError);
            }
        }
    }


    // Get average position of boundary face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    vectorField avgBoundary(pointFaces.size(), vector::zero);
    labelList nBoundary(pointFaces.size(), 0);

    forAll(pointFaces, patchPointI)
    {
        const labelList& pFaces = pointFaces[patchPointI];

        forAll(pFaces, pfI)
        {
            label faceI = pFaces[pfI];

            if (isMasterFace.get(pp.addressing()[faceI]) == 1)
            {
                avgBoundary[patchPointI] += pp[faceI].centre(points);
                nBoundary[patchPointI]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        avgBoundary,
        plusEqOp<point>(),  // combine op
        vector::zero,       // null value
        false               // no separation
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        nBoundary,
        plusEqOp<label>(),  // combine op
        0,                  // null value
        false               // no separation
    );

    forAll(avgBoundary, i)
    {
        avgBoundary[i] /= nBoundary[i];
    }


    // Get average position of internal face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    vectorField avgInternal;
    labelList nInternal;
    {
        vectorField globalSum(mesh.nPoints(), vector::zero);
        labelList globalNum(mesh.nPoints(), 0);

        // Note: no use of pointFaces
        const faceList& faces = mesh.faces();

        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            const face& f = faces[faceI];
            const point& fc = mesh.faceCentres()[faceI];

            forAll(f, fp)
            {
                globalSum[f[fp]] += fc;
                globalNum[f[fp]]++;
            }
        }

        // Count coupled faces as internal ones (but only once)
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        forAll(patches, patchI)
        {
            if (Pstream::parRun() && isA<processorPolyPatch>(patches[patchI]))
            {
                const processorPolyPatch& pp =
                    refCast<const processorPolyPatch>(patches[patchI]);

                if (pp.myProcNo() < pp.neighbProcNo())
                {
                    const vectorField::subField faceCentres = pp.faceCentres();

                    forAll(pp, i)
                    {
                        const face& f = pp[i];
                        const point& fc = faceCentres[i];

                        forAll(f, fp)
                        {
                            globalSum[f[fp]] += fc;
                            globalNum[f[fp]]++;
                        }
                    }
                }
            }
            else if (isA<cyclicPolyPatch>(patches[patchI]))
            {
                const cyclicPolyPatch& pp =
                    refCast<const cyclicPolyPatch>(patches[patchI]);

                const vectorField::subField faceCentres = pp.faceCentres();

                for (label i = 0; i < pp.size()/2; i++)
                {
                    const face& f = pp[i];
                    const point& fc = faceCentres[i];

                    forAll(f, fp)
                    {
                        globalSum[f[fp]] += fc;
                        globalNum[f[fp]]++;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            globalSum,
            plusEqOp<vector>(), // combine op
            vector::zero,       // null value
            false               // no separation
        );
        syncTools::syncPointList
        (
            mesh,
            globalNum,
            plusEqOp<label>(),  // combine op
            0,                  // null value
            false               // no separation
        );

        avgInternal.setSize(meshPoints.size());
        nInternal.setSize(meshPoints.size());

        forAll(avgInternal, patchPointI)
        {
            label meshPointI = meshPoints[patchPointI];

            nInternal[patchPointI] = globalNum[meshPointI];

            if (nInternal[patchPointI] == 0)
            {
                avgInternal[patchPointI] = globalSum[meshPointI];
            }
            else
            {
                avgInternal[patchPointI] =
                    globalSum[meshPointI]
                  / nInternal[patchPointI];
            }
        }
    }


    // Precalculate any cell using mesh point (replacement of pointCells()[])
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList anyCell(mesh.nPoints(), -1);
    forAll(mesh.faceNeighbour(), faceI)
    {
        label own = mesh.faceOwner()[faceI];
        const face& f = mesh.faces()[faceI];

        forAll(f, fp)
        {
            anyCell[f[fp]] = own;
        }
    }
    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        label own = mesh.faceOwner()[faceI];

        const face& f = mesh.faces()[faceI];

        forAll(f, fp)
        {
            anyCell[f[fp]] = own;
        }
    }


    // Displacement to calculate.
    pointField patchDisp(meshPoints.size(), vector::zero);

    forAll(pointFaces, i)
    {
        label meshPointI = meshPoints[i];
        const point& currentPos = pp.points()[meshPointI];

        // Now we have the two average points: avgBoundary and avgInternal
        // and how many boundary/internal faces connect to the point
        // (nBoundary, nInternal)
        // Do some blending between the two.
        // Note: the following section has some reasoning behind it but the
        // blending factors can be experimented with.

        point newPos;

        if (nonManifoldPoint.get(i) == 0u)
        {
            // Points that are manifold. Weight the internal and boundary
            // by their number of faces and blend with
            scalar internalBlend = 0.1;
            scalar blend = 0.1;

            point avgPos =
                (
                   internalBlend*nInternal[i]*avgInternal[i]
                  +(1-internalBlend)*nBoundary[i]*avgBoundary[i]
                )
              / (internalBlend*nInternal[i]+(1-internalBlend)*nBoundary[i]);

            newPos = (1-blend)*avgPos + blend*currentPos;
        }
        else if (nInternal[i] == 0)
        {
            // Non-manifold without internal faces. Use any connected cell
            // as internal point instead. Use precalculated any cell to avoid
            // e.g. pointCells()[meshPointI][0]

            const point& cc = mesh.cellCentres()[anyCell[meshPointI]];

            scalar cellCBlend = 0.8;
            scalar blend = 0.1;

            point avgPos = (1-cellCBlend)*avgBoundary[i] + cellCBlend*cc;

            newPos = (1-blend)*avgPos + blend*currentPos;
        }
        else
        {
            // Non-manifold point with internal faces connected to them
            scalar internalBlend = 0.9;
            scalar blend = 0.1;

            point avgPos =
                internalBlend*avgInternal[i]
              + (1-internalBlend)*avgBoundary[i];

            newPos = (1-blend)*avgPos + blend*currentPos;
        }

        patchDisp[i] = newPos - currentPos;
    }

    return patchDisp;
}


Foam::tmp<Foam::scalarField> Foam::autoSnapDriver::edgePatchDist
(
    const pointMesh& pMesh,
    const indirectPrimitivePatch& pp
)
{
    const polyMesh& mesh = pMesh();

    // Set initial changed points to all the patch points
    List<pointEdgePoint> wallInfo(pp.nPoints());

    forAll(pp.localPoints(), ppI)
    {
        wallInfo[ppI] = pointEdgePoint(pp.localPoints()[ppI], 0.0);
    }

    // Current info on points
    List<pointEdgePoint> allPointInfo(mesh.nPoints());

    // Current info on edges
    List<pointEdgePoint> allEdgeInfo(mesh.nEdges());

    PointEdgeWave<pointEdgePoint> wallCalc
    (
        mesh,
        pp.meshPoints(),
        wallInfo,

        allPointInfo,
        allEdgeInfo,
        mesh.globalData().nTotalPoints()  // max iterations
    );

    // Copy edge values into scalarField
    tmp<scalarField> tedgeDist(new scalarField(mesh.nEdges()));
    scalarField& edgeDist = tedgeDist();

    forAll(allEdgeInfo, edgeI)
    {
        edgeDist[edgeI] = Foam::sqrt(allEdgeInfo[edgeI].distSqr());
    }


    //{
    //    // For debugging: dump to file
    //    pointScalarField pointDist
    //    (
    //        IOobject
    //        (
    //            "pointDist",
    //            meshRefiner_.timeName(),
    //            mesh.DB(),
    //            IOobject::NO_READ,
    //            IOobject::AUTO_WRITE
    //        ),
    //        pMesh,
    //        dimensionedScalar("pointDist", dimless, 0.0)
    //    );
    //
    //    forAll(allEdgeInfo, edgeI)
    //    {
    //        scalar d = Foam::sqrt(allEdgeInfo[edgeI].distSqr());
    //
    //        const edge& e = mesh.edges()[edgeI];
    //
    //        pointDist[e[0]] += d;
    //        pointDist[e[1]] += d;
    //    }
    //    forAll(pointDist, pointI)
    //    {
    //        pointDist[pointI] /= mesh.pointEdges()[pointI].size();
    //    }
    //    Info<< "Writing patch distance to " << pointDist.name()
    //        << " at time " << meshRefiner_.timeName() << endl;
    //
    //    pointDist.write();
    //}

    return tedgeDist;
}


void Foam::autoSnapDriver::dumpMove
(
    const fileName& fName,
    const pointField& meshPts,
    const pointField& surfPts
)
{
    // Dump direction of growth into file
    Pout<< nl << "Dumping move direction to " << fName << nl
        << "View this Lightwave-OBJ file with e.g. javaview" << nl
        << endl;

    OFstream nearestStream(fName);

    label vertI = 0;

    forAll(meshPts, ptI)
    {
        meshTools::writeOBJ(nearestStream, meshPts[ptI]);
        vertI++;

        meshTools::writeOBJ(nearestStream, surfPts[ptI]);
        vertI++;

        nearestStream<< "l " << vertI-1 << ' ' << vertI << nl;
    }
}


// Check whether all displacement vectors point outwards of patch. Return true
// if so.
bool Foam::autoSnapDriver::outwardsDisplacement
(
    const indirectPrimitivePatch& pp,
    const vectorField& patchDisp
)
{
    const vectorField& faceNormals = pp.faceNormals();
    const labelListList& pointFaces = pp.pointFaces();

    forAll(pointFaces, pointI)
    {
        const labelList& pFaces = pointFaces[pointI];

        vector disp(patchDisp[pointI]);

        scalar magDisp = mag(disp);

        if (magDisp > SMALL)
        {
            disp /= magDisp;

            bool outwards = meshTools::visNormal(disp, faceNormals, pFaces);

            if (!outwards)
            {
                Warning<< "Displacement " << patchDisp[pointI]
                    << " at mesh point " << pp.meshPoints()[pointI]
                    << " coord " << pp.points()[pp.meshPoints()[pointI]]
                    << " points through the surrounding patch faces" << endl;
                return false;
            }
        }
        else
        {
            //? Displacement small but in wrong direction. Would probably be ok.
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoSnapDriver::autoSnapDriver
(
    meshRefinement& meshRefiner,
    const labelList& globalToPatch
)
:
    meshRefiner_(meshRefiner),
    globalToPatch_(globalToPatch)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::mapPolyMesh> Foam::autoSnapDriver::createZoneBaffles
(
    List<labelPair>& baffles
)
{
    labelList zonedSurfaces = meshRefiner_.surfaces().getNamedSurfaces();

    autoPtr<mapPolyMesh> map;

    // No need to sync; all processors will have all same zonedSurfaces.
    if (zonedSurfaces.size())
    {
        fvMesh& mesh = meshRefiner_.mesh();

        // Split internal faces on interface surfaces
        Info<< "Converting zoned faces into baffles ..." << endl;

        // Get faces (internal only) to be baffled. Map from face to patch
        // label.
        Map<label> faceToPatch(getZoneBafflePatches(false));

        label nZoneFaces = returnReduce(faceToPatch.size(), sumOp<label>());
        if (nZoneFaces > 0)
        {
            // Convert into labelLists
            labelList ownPatch(mesh.nFaces(), -1);
            forAllConstIter(Map<label>, faceToPatch, iter)
            {
                ownPatch[iter.key()] = iter();
            }

            // Create baffles. both sides same patch.
            map = meshRefiner_.createBaffles(ownPatch, ownPatch);

            // Get pairs of faces created.
            // Just loop over faceMap and store baffle if we encounter a slave
            // face.

            baffles.setSize(faceToPatch.size());
            label baffleI = 0;

            const labelList& faceMap = map().faceMap();
            const labelList& reverseFaceMap = map().reverseFaceMap();

            forAll(faceMap, faceI)
            {
                label oldFaceI = faceMap[faceI];

                // Does face originate from face-to-patch
                Map<label>::const_iterator iter = faceToPatch.find(oldFaceI);

                if (iter != faceToPatch.end())
                {
                    label masterFaceI = reverseFaceMap[oldFaceI];
                    if (faceI != masterFaceI)
                    {
                        baffles[baffleI++] = labelPair(masterFaceI, faceI);
                    }
                }
            }

            if (baffleI != faceToPatch.size())
            {
                FatalErrorIn("autoSnapDriver::createZoneBaffles(..)")
                    << "Had " << faceToPatch.size() << " patches to create "
                    << " but encountered " << baffleI
                    << " slave faces originating from patcheable faces."
                    << abort(FatalError);
            }

            if (debug)
            {
                const_cast<Time&>(mesh.time())++;
                Pout<< "Writing baffled mesh to time "
                    << meshRefiner_.timeName() << endl;
                mesh.write();
            }
        }
        Info<< "Created " << nZoneFaces << " baffles in = "
            << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }
    return map;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::autoSnapDriver::mergeZoneBaffles
(
    const List<labelPair>& baffles
)
{
    labelList zonedSurfaces = meshRefiner_.surfaces().getNamedSurfaces();

    autoPtr<mapPolyMesh> map;

    // No need to sync; all processors will have all same zonedSurfaces.
    label nBaffles = returnReduce(baffles.size(), sumOp<label>());
    if (zonedSurfaces.size() && nBaffles > 0)
    {
        // Merge any baffles
        Info<< "Converting " << nBaffles << " baffles back into zoned faces ..."
            << endl;

        map = meshRefiner_.mergeBaffles(baffles);

        Info<< "Converted baffles in = "
            << meshRefiner_.mesh().time().cpuTimeIncrement()
            << " s\n" << nl << endl;
    }

    return map;
}


Foam::scalarField Foam::autoSnapDriver::calcSnapDistance
(
    const snapParameters& snapParams,
    const indirectPrimitivePatch& pp
) const
{
    const edgeList& edges = pp.edges();
    const labelListList& pointEdges = pp.pointEdges();
    const pointField& localPoints = pp.localPoints();
    const fvMesh& mesh = meshRefiner_.mesh();

    scalarField maxEdgeLen(localPoints.size(), -GREAT);

    forAll(pointEdges, pointI)
    {
        const labelList& pEdges = pointEdges[pointI];

        forAll(pEdges, pEdgeI)
        {
            const edge& e = edges[pEdges[pEdgeI]];

            scalar len = e.mag(localPoints);

            maxEdgeLen[pointI] = max(maxEdgeLen[pointI], len);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        maxEdgeLen,
        maxEqOp<scalar>(),  // combine op
        -GREAT,             // null value
        false               // no separation
    );

    return snapParams.snapTol()*maxEdgeLen;
}


void Foam::autoSnapDriver::preSmoothPatch
(
    const snapParameters& snapParams,
    const label nInitErrors,
    const List<labelPair>& baffles,
    motionSmoother& meshMover
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    labelList checkFaces;

    Info<< "Smoothing patch points ..." << endl;
    for
    (
        label smoothIter = 0;
        smoothIter < snapParams.nSmoothPatch();
        smoothIter++
    )
    {
        Info<< "Smoothing iteration " << smoothIter << endl;
        checkFaces.setSize(mesh.nFaces());
        forAll(checkFaces, faceI)
        {
            checkFaces[faceI] = faceI;
        }

        pointField patchDisp(smoothPatchDisplacement(meshMover, baffles));

        // The current mesh is the starting mesh to smooth from.
        meshMover.setDisplacement(patchDisp);
        meshMover.correct();

        scalar oldErrorReduction = -1;

        for (label snapIter = 0; snapIter < 2*snapParams.nSnap(); snapIter++)
        {
            Info<< nl << "Scaling iteration " << snapIter << endl;

            if (snapIter == snapParams.nSnap())
            {
                Info<< "Displacement scaling for error reduction set to 0."
                    << endl;
                oldErrorReduction = meshMover.setErrorReduction(0.0);
            }

            // Try to adapt mesh to obtain displacement by smoothly
            // decreasing displacement at error locations.
            if (meshMover.scaleMesh(checkFaces, baffles, true, nInitErrors))
            {
                Info<< "Successfully moved mesh" << endl;
                break;
            }
        }

        if (oldErrorReduction >= 0)
        {
            meshMover.setErrorReduction(oldErrorReduction);
        }
        Info<< endl;
    }


    // The current mesh is the starting mesh to smooth from.
    meshMover.correct();

    if (debug)
    {
        const_cast<Time&>(mesh.time())++;
        Pout<< "Writing patch smoothed mesh to time " << meshRefiner_.timeName()
            << endl;

        mesh.write();
    }

    Info<< "Patch points smoothed in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
}


// Get (pp-local) indices of points that are both on zone and on patched surface
Foam::labelList Foam::autoSnapDriver::getZoneSurfacePoints
(
    const indirectPrimitivePatch& pp,
    const word& zoneName
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    label zoneI = mesh.faceZones().findZoneID(zoneName);

    if (zoneI == -1)
    {
        FatalErrorIn
        (
            "autoSnapDriver::getZoneSurfacePoints"
            "(const indirectPrimitivePatch&, const word&)"
        )   << "Cannot find zone " << zoneName
            << exit(FatalError);
    }

    const faceZone& fZone = mesh.faceZones()[zoneI];


    // Could use PrimitivePatch & localFaces to extract points but might just
    // as well do it ourselves.

    boolList pointOnZone(pp.nPoints(), false);

    forAll(fZone, i)
    {
        const face& f = mesh.faces()[fZone[i]];

        forAll(f, fp)
        {
            label meshPointI = f[fp];

            Map<label>::const_iterator iter =
                pp.meshPointMap().find(meshPointI);

            if (iter != pp.meshPointMap().end())
            {
                label pointI = iter();
                pointOnZone[pointI] = true;
            }
        }
    }

    return findIndices(pointOnZone, true);
}


Foam::vectorField Foam::autoSnapDriver::calcNearestSurface
(
    const scalarField& snapDist,
    motionSmoother& meshMover
) const
{
    Info<< "Calculating patchDisplacement as distance to nearest surface"
        << " point ..." << endl;

    const indirectPrimitivePatch& pp = meshMover.patch();
    const pointField& localPoints = pp.localPoints();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();
    const fvMesh& mesh = meshRefiner_.mesh();

    // Displacement per patch point
    vectorField patchDisp(localPoints.size(), vector::zero);

    if (returnReduce(localPoints.size(), sumOp<label>()) > 0)
    {
        // Current surface snapped to
        labelList snapSurf(localPoints.size(), -1);

        // Divide surfaces into zoned and unzoned
        labelList zonedSurfaces =
            meshRefiner_.surfaces().getNamedSurfaces();
        labelList unzonedSurfaces =
            meshRefiner_.surfaces().getUnnamedSurfaces();


        // 1. All points to non-interface surfaces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        {
            List<pointIndexHit> hitInfo;
            labelList hitSurface;
            surfaces.findNearest
            (
                unzonedSurfaces,
                localPoints,
                sqr(4*snapDist),        // sqr of attract distance
                hitSurface,
                hitInfo
            );

            forAll(hitInfo, pointI)
            {
                if (hitInfo[pointI].hit())
                {
                    patchDisp[pointI] =
                        hitInfo[pointI].hitPoint()
                      - localPoints[pointI];

                    snapSurf[pointI] = hitSurface[pointI];
                }
            }
        }



        // 2. All points on zones to their respective surface
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Surfaces with zone information
        const wordList& faceZoneNames = surfaces.faceZoneNames();

        // Current best snap distance
        scalarField minSnapDist(snapDist);

        forAll(zonedSurfaces, i)
        {
            label zoneSurfI = zonedSurfaces[i];

            const labelList surfacesToTest(1, zoneSurfI);

            // Get indices of points both on faceZone and on pp.
            labelList zonePointIndices
            (
                getZoneSurfacePoints
                (
                    pp,
                    faceZoneNames[zoneSurfI]
                )
            );

            // Find nearest for points both on faceZone and pp.
            List<pointIndexHit> hitInfo;
            labelList hitSurface;
            surfaces.findNearest
            (
                labelList(1, zoneSurfI),
                pointField(localPoints, zonePointIndices),
                sqr(4*scalarField(minSnapDist, zonePointIndices)),
                hitSurface,
                hitInfo
            );

            forAll(hitInfo, i)
            {
                label pointI = zonePointIndices[i];

                if (hitInfo[i].hit())
                {
                    patchDisp[pointI] =
                        hitInfo[i].hitPoint()
                      - localPoints[pointI];

                    minSnapDist[pointI] = min
                    (
                        minSnapDist[pointI],
                        mag(patchDisp[pointI])
                    );

                    snapSurf[pointI] = zoneSurfI;
                }
            }
        }

        // Check if all points are being snapped
        forAll(snapSurf, pointI)
        {
            if (snapSurf[pointI] == -1)
            {
                WarningIn("autoSnapDriver::calcNearestSurface(..)")
                    << "For point:" << pointI
                    << " coordinate:" << localPoints[pointI]
                    << " did not find any surface within:"
                    << minSnapDist[pointI]
                    << " meter." << endl;
            }
        }

        {
            scalarField magDisp(mag(patchDisp));

            Info<< "Wanted displacement : average:"
                << gSum(magDisp)/returnReduce(patchDisp.size(), sumOp<label>())
                << " min:" << gMin(magDisp)
                << " max:" << gMax(magDisp) << endl;
        }
    }

    Info<< "Calculated surface displacement in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;


    // Limit amount of movement.
    forAll(patchDisp, patchPointI)
    {
        scalar magDisp = mag(patchDisp[patchPointI]);

        if (magDisp > snapDist[patchPointI])
        {
            patchDisp[patchPointI] *= snapDist[patchPointI] / magDisp;

            Pout<< "Limiting displacement for " << patchPointI
                << " from " << magDisp << " to " << snapDist[patchPointI]
                << endl;
        }
    }

    // Points on zones in one domain but only present as point on other
    // will not do condition 2 on all. Sync explicitly.
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        patchDisp,
        minMagEqOp(),                   // combine op
        vector(GREAT, GREAT, GREAT),    // null value
        false                           // no separation
    );


    // Check for displacement being outwards.
    outwardsDisplacement(pp, patchDisp);

    // Set initial distribution of displacement field (on patches) from
    // patchDisp and make displacement consistent with b.c. on displacement
    // pointVectorField.
    meshMover.setDisplacement(patchDisp);

    if (debug)
    {
        dumpMove
        (
            mesh.time().path()/"patchDisplacement.obj",
            pp.localPoints(),
            pp.localPoints() + patchDisp
        );
    }

    return patchDisp;
}


void Foam::autoSnapDriver::smoothDisplacement
(
    const snapParameters& snapParams,
    motionSmoother& meshMover
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const pointMesh& pMesh = meshMover.pMesh();
    const indirectPrimitivePatch& pp = meshMover.patch();

    Info<< "Smoothing displacement ..." << endl;

    // Set edge diffusivity as inverse of distance to patch
    scalarField edgeGamma(1.0/(edgePatchDist(pMesh, pp) + SMALL));
    //scalarField edgeGamma(mesh.nEdges(), 1.0);
    //scalarField edgeGamma(wallGamma(mesh, pp, 10, 1));

    // Get displacement field
    pointVectorField& disp = meshMover.displacement();

    for (label iter = 0; iter < snapParams.nSmoothDispl(); iter++)
    {
        if ((iter % 10) == 0)
        {
            Info<< "Iteration " << iter << endl;
        }
        pointVectorField oldDisp(disp);

        meshMover.smooth(oldDisp, edgeGamma, false, disp);
    }
    Info<< "Displacement smoothed in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;

    if (debug)
    {
        const_cast<Time&>(mesh.time())++;
        Pout<< "Writing smoothed mesh to time " << meshRefiner_.timeName()
            << endl;

        // Moving mesh creates meshPhi. Can be cleared out by a mesh.clearOut
        // but this will also delete all pointMesh but not pointFields which
        // gives an illegal situation.

        mesh.write();

        Pout<< "Writing displacement field ..." << endl;
        disp.write();
        tmp<pointScalarField> magDisp(mag(disp));
        magDisp().write();

        Pout<< "Writing actual patch displacement ..." << endl;
        vectorField actualPatchDisp(disp, pp.meshPoints());
        dumpMove
        (
            mesh.time().path()/"actualPatchDisplacement.obj",
            pp.localPoints(),
            pp.localPoints() + actualPatchDisp
        );
    }
}


void Foam::autoSnapDriver::scaleMesh
(
    const snapParameters& snapParams,
    const label nInitErrors,
    const List<labelPair>& baffles,
    motionSmoother& meshMover
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    // Relax displacement until correct mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    labelList checkFaces(identity(mesh.nFaces()));

    scalar oldErrorReduction = -1;

    Info<< "Moving mesh ..." << endl;
    for (label iter = 0; iter < 2*snapParams.nSnap(); iter++)
    {
        Info<< nl << "Iteration " << iter << endl;

        if (iter == snapParams.nSnap())
        {
            Info<< "Displacement scaling for error reduction set to 0." << endl;
            oldErrorReduction = meshMover.setErrorReduction(0.0);
        }

        if (meshMover.scaleMesh(checkFaces, baffles, true, nInitErrors))
        {
            Info<< "Successfully moved mesh" << endl;

            break;
        }
        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
            Pout<< "Writing scaled mesh to time " << meshRefiner_.timeName()
                << endl;
            mesh.write();

            Pout<< "Writing displacement field ..." << endl;
            meshMover.displacement().write();
            tmp<pointScalarField> magDisp(mag(meshMover.displacement()));
            magDisp().write();
        }
    }

    if (oldErrorReduction >= 0)
    {
        meshMover.setErrorReduction(oldErrorReduction);
    }
    Info<< "Moved mesh in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
}


// After snapping: correct patching according to nearest surface.
// Code is very similar to calcNearestSurface.
// - calculate face-wise snap distance as max of point-wise
// - calculate face-wise nearest surface point
// - repatch face according to patch for surface point.
Foam::autoPtr<Foam::mapPolyMesh> Foam::autoSnapDriver::repatchToSurface
(
    const snapParameters& snapParams,
    const labelList& adaptPatchIDs
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    Info<< "Repatching faces according to nearest surface ..." << endl;

    // Get the labels of added patches.
    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            adaptPatchIDs
        )
    );
    indirectPrimitivePatch& pp = ppPtr();

    // Divide surfaces into zoned and unzoned
    labelList zonedSurfaces = meshRefiner_.surfaces().getNamedSurfaces();
    labelList unzonedSurfaces = meshRefiner_.surfaces().getUnnamedSurfaces();


    // Faces that do not move
    PackedBoolList isZonedFace(mesh.nFaces(), 0);
    {
        // 1. All faces on zoned surfaces
        const wordList& faceZoneNames = surfaces.faceZoneNames();
        const faceZoneMesh& fZones = mesh.faceZones();

        forAll(zonedSurfaces, i)
        {
            label zoneSurfI = zonedSurfaces[i];

            label zoneI = fZones.findZoneID(faceZoneNames[zoneSurfI]);

            const faceZone& fZone = fZones[zoneI];

            forAll(fZone, i)
            {
                isZonedFace.set(fZone[i], 1);
            }
        }
    }


    // Determine per pp face which patch it should be in
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Patch that face should be in
    labelList closestPatch(pp.size(), -1);
    {
        // face snap distance as max of point snap distance
        scalarField faceSnapDist(pp.size(), -GREAT);
        {
            // Distance to attract to nearest feature on surface
            const scalarField snapDist(calcSnapDistance(snapParams, pp));

            const faceList& localFaces = pp.localFaces();

            forAll(localFaces, faceI)
            {
                const face& f = localFaces[faceI];

                forAll(f, fp)
                {
                    faceSnapDist[faceI] = max
                    (
                        faceSnapDist[faceI],
                        snapDist[f[fp]]
                    );
                }
            }
        }

        pointField localFaceCentres(mesh.faceCentres(), pp.addressing());

        // Get nearest surface and region
        labelList hitSurface;
        labelList hitRegion;
        surfaces.findNearestRegion
        (
            unzonedSurfaces,
            localFaceCentres,
            sqr(4*faceSnapDist),    // sqr of attract distance
            hitSurface,
            hitRegion
        );

        // Get patch
        forAll(pp, i)
        {
            label faceI = pp.addressing()[i];

            if (hitSurface[i] != -1 && (isZonedFace.get(faceI) == 0))
            {
                closestPatch[i] = globalToPatch_
                [
                    surfaces.globalRegion
                    (
                        hitSurface[i],
                        hitRegion[i]
                    )
                ];
            }
        }
    }


    // Change those faces for which there is a different closest patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList ownPatch(mesh.nFaces(), -1);
    labelList neiPatch(mesh.nFaces(), -1);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        forAll(pp, i)
        {
            ownPatch[pp.start()+i] = patchI;
            neiPatch[pp.start()+i] = patchI;
        }
    }

    label nChanged = 0;
    forAll(closestPatch, i)
    {
        label faceI = pp.addressing()[i];

        if (closestPatch[i] != -1 && closestPatch[i] != ownPatch[faceI])
        {
            ownPatch[faceI] = closestPatch[i];
            neiPatch[faceI] = closestPatch[i];
            nChanged++;
        }
    }

    Info<< "Repatched " << returnReduce(nChanged, sumOp<label>())
        << " faces in = " << mesh.time().cpuTimeIncrement() << " s\n" << nl
        << endl;

    return meshRefiner_.createBaffles(ownPatch, neiPatch);
}


void Foam::autoSnapDriver::doSnap
(
    const dictionary& snapDict,
    const dictionary& motionDict,
    const snapParameters& snapParams
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl
        << "Morphing phase" << nl
        << "--------------" << nl
        << endl;

    // Get the labels of added patches.
    labelList adaptPatchIDs(meshRefiner_.meshedPatches());

    // Create baffles (pairs of faces that share the same points)
    // Baffles stored as owner and neighbour face that have been created.
    List<labelPair> baffles;
    createZoneBaffles(baffles);

    {
        autoPtr<indirectPrimitivePatch> ppPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                adaptPatchIDs
            )
        );
        indirectPrimitivePatch& pp = ppPtr();

        // Distance to attract to nearest feature on surface
        const scalarField snapDist(calcSnapDistance(snapParams, pp));


        // Construct iterative mesh mover.
        Info<< "Constructing mesh displacer ..." << endl;
        Info<< "Using mesh parameters " << motionDict << nl << endl;

        const pointMesh& pMesh = pointMesh::New(mesh);

        motionSmoother meshMover
        (
            mesh,
            pp,
            adaptPatchIDs,
            meshRefinement::makeDisplacementField(pMesh, adaptPatchIDs),
            motionDict
        );


        // Check initial mesh
        Info<< "Checking initial mesh ..." << endl;
        labelHashSet wrongFaces(mesh.nFaces()/100);
        motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces);
        const label nInitErrors = returnReduce
        (
            wrongFaces.size(),
            sumOp<label>()
        );

        Info<< "Detected " << nInitErrors << " illegal faces"
            << " (concave, zero area or negative cell pyramid volume)"
            << endl;


        Info<< "Checked initial mesh in = "
            << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;

        // Pre-smooth patch vertices (so before determining nearest)
        preSmoothPatch(snapParams, nInitErrors, baffles, meshMover);

        // Calculate displacement at every patch point. Insert into
        // meshMover.
        calcNearestSurface(snapDist, meshMover);

        //// Get smoothly varying internal displacement field.
        //- 2009-12-16 : was not found to be beneficial. Keeping internal
        // fields fixed slightly increases skewness (on boundary)
        // but lowers non-orthogonality quite a bit (e.g. 65->59 degrees).
        // Maybe if better smoother?
        //smoothDisplacement(snapParams, meshMover);

        // Apply internal displacement to mesh.
        scaleMesh(snapParams, nInitErrors, baffles, meshMover);
    }

    // Merge any introduced baffles.
    mergeZoneBaffles(baffles);

    // Repatch faces according to nearest.
    repatchToSurface(snapParams, adaptPatchIDs);
}


// ************************************************************************* //
