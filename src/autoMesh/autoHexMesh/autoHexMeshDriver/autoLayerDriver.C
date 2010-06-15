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
    All to do with adding cell layers

\*----------------------------------------------------------------------------*/

#include "autoLayerDriver.H"
#include "fvMesh.H"
#include "Time.H"
#include "meshRefinement.H"
#include "removePoints.H"
#include "pointFields.H"
#include "motionSmoother.H"
#include "mathematicalConstants.H"
#include "pointSet.H"
#include "faceSet.H"
#include "cellSet.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "addPatchCellLayer.H"
#include "mapDistributePolyMesh.H"
#include "OFstream.H"
#include "layerParameters.H"
#include "combineFaces.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(autoLayerDriver, 0);

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::autoLayerDriver::mergePatchFacesUndo
(
    const scalar minCos,
    const scalar concaveCos,
    const dictionary& motionDict
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // Patch face merging engine
    combineFaces faceCombiner(mesh, true);

    // Pick up all candidate cells on boundary
    labelHashSet boundaryCells(mesh.nFaces()-mesh.nInternalFaces());

    {
        labelList patchIDs(meshRefiner_.meshedPatches());

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        forAll(patchIDs, i)
        {
            label patchI = patchIDs[i];

            const polyPatch& patch = patches[patchI];

            if (!patch.coupled())
            {
                forAll(patch, i)
                {
                    boundaryCells.insert(mesh.faceOwner()[patch.start()+i]);
                }
            }
        }
    }

    // Get all sets of faces that can be merged
    labelListList allFaceSets
    (
        faceCombiner.getMergeSets
        (
            minCos,
            concaveCos,
            boundaryCells
        )
    );

    label nFaceSets = returnReduce(allFaceSets.size(), sumOp<label>());

    Info<< "Merging " << nFaceSets << " sets of faces." << nl << endl;

    if (nFaceSets > 0)
    {
        if (debug)
        {
            faceSet allSets(mesh, "allFaceSets", allFaceSets.size());
            forAll(allFaceSets, setI)
            {
                forAll(allFaceSets[setI], i)
                {
                    allSets.insert(allFaceSets[setI][i]);
                }
            }
            Pout<< "Writing all faces to be merged to set "
                << allSets.objectPath() << endl;
            allSets.write();
        }


        // Topology changes container
        polyTopoChange meshMod(mesh);

        // Merge all faces of a set into the first face of the set.
        faceCombiner.setRefinement(allFaceSets, meshMod);

        // Experimental: store data for all the points that have been deleted
        meshRefiner_.storeData
        (
            faceCombiner.savedPointLabels(),    // points to store
            labelList(0),                       // faces to store
            labelList(0)                        // cells to store
        );

        // Change the mesh (no inflation)
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

        // Update fields
        mesh.updateMesh(map);

        // Move mesh (since morphing does not do this)
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh.clearOut();
        }

        if (meshRefiner_.overwrite())
        {
            mesh.setInstance(meshRefiner_.oldInstance());
        }

        faceCombiner.updateMesh(map);

        meshRefiner_.updateMesh(map, labelList(0));


        for (label iteration = 0; iteration < 100; iteration++)
        {
            Info<< nl
                << "Undo iteration " << iteration << nl
                << "----------------" << endl;


            // Check mesh for errors
            // ~~~~~~~~~~~~~~~~~~~~~

            faceSet errorFaces
            (
                mesh,
                "errorFaces",
                mesh.nFaces()-mesh.nInternalFaces()
            );
            bool hasErrors = motionSmoother::checkMesh
            (
                false,  // report
                mesh,
                motionDict,
                errorFaces
            );
            //if (checkEdgeConnectivity)
            //{
            //    Info<< "Checking edge-face connectivity (duplicate faces"
            //        << " or non-consecutive shared vertices)" << endl;
            //
            //    label nOldSize = errorFaces.size();
            //
            //    hasErrors =
            //        mesh.checkFaceFaces
            //        (
            //            false,
            //            &errorFaces
            //        )
            //     || hasErrors;
            //
            //    Info<< "Detected additional "
            //        << returnReduce(errorFaces.size()-nOldSize, sumOp<label>())
            //        << " faces with illegal face-face connectivity"
            //        << endl;
            //}

            if (!hasErrors)
            {
                break;
            }


            if (debug)
            {
                Pout<< "Writing all faces in error to faceSet "
                    << errorFaces.objectPath() << nl << endl;
                errorFaces.write();
            }


            // Check any master cells for using any of the error faces
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            DynamicList<label> mastersToRestore(allFaceSets.size());

            forAll(allFaceSets, setI)
            {
                label masterFaceI = faceCombiner.masterFace()[setI];

                if (masterFaceI != -1)
                {
                    label masterCellII = mesh.faceOwner()[masterFaceI];

                    const cell& cFaces = mesh.cells()[masterCellII];

                    forAll(cFaces, i)
                    {
                        if (errorFaces.found(cFaces[i]))
                        {
                            mastersToRestore.append(masterFaceI);
                            break;
                        }
                    }
                }
            }
            mastersToRestore.shrink();

            label nRestore = returnReduce
            (
                mastersToRestore.size(),
                sumOp<label>()
            );

            Info<< "Masters that need to be restored:"
                << nRestore << endl;

            if (debug)
            {
                faceSet restoreSet(mesh, "mastersToRestore", mastersToRestore);
                Pout<< "Writing all " << mastersToRestore.size()
                    << " masterfaces to be restored to set "
                    << restoreSet.objectPath() << endl;
                restoreSet.write();
            }


            if (nRestore == 0)
            {
                break;
            }


            // Undo
            // ~~~~

            // Topology changes container
            polyTopoChange meshMod(mesh);

            // Merge all faces of a set into the first face of the set.
            // Experimental:mark all points/faces/cells that have been restored.
            Map<label> restoredPoints(0);
            Map<label> restoredFaces(0);
            Map<label> restoredCells(0);

            faceCombiner.setUnrefinement
            (
                mastersToRestore,
                meshMod,
                restoredPoints,
                restoredFaces,
                restoredCells
            );

            // Change the mesh (no inflation)
            autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

            // Update fields
            mesh.updateMesh(map);

            // Move mesh (since morphing does not do this)
            if (map().hasMotionPoints())
            {
                mesh.movePoints(map().preMotionPoints());
            }
            else
            {
                // Delete mesh volumes.
                mesh.clearOut();
            }

            if (meshRefiner_.overwrite())
            {
                mesh.setInstance(meshRefiner_.oldInstance());
            }

            faceCombiner.updateMesh(map);

            // Renumber restore maps
            inplaceMapKey(map().reversePointMap(), restoredPoints);
            inplaceMapKey(map().reverseFaceMap(), restoredFaces);
            inplaceMapKey(map().reverseCellMap(), restoredCells);

            // Experimental:restore all points/face/cells in maps
            meshRefiner_.updateMesh
            (
                map,
                labelList(0),       // changedFaces
                restoredPoints,
                restoredFaces,
                restoredCells
            );

            Info<< endl;
        }

        if (debug)
        {
            Pout<< "Writing merged-faces mesh to time "
                << meshRefiner_.timeName() << nl << endl;
            mesh.write();
        }
    }
    else
    {
        Info<< "No faces merged ..." << endl;
    }

    return nFaceSets;
}


// Remove points. pointCanBeDeleted is parallel synchronised.
Foam::autoPtr<Foam::mapPolyMesh> Foam::autoLayerDriver::doRemovePoints
(
    removePoints& pointRemover,
    const boolList& pointCanBeDeleted
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // Topology changes container
    polyTopoChange meshMod(mesh);

    pointRemover.setRefinement(pointCanBeDeleted, meshMod);

    // Change the mesh (no inflation)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh.clearOut();
    }

    if (meshRefiner_.overwrite())
    {
        mesh.setInstance(meshRefiner_.oldInstance());
    }

    pointRemover.updateMesh(map);
    meshRefiner_.updateMesh(map, labelList(0));

    return map;
}


// Restore faces (which contain removed points)
Foam::autoPtr<Foam::mapPolyMesh> Foam::autoLayerDriver::doRestorePoints
(
    removePoints& pointRemover,
    const labelList& facesToRestore
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // Topology changes container
    polyTopoChange meshMod(mesh);

    // Determine sets of points and faces to restore
    labelList localFaces, localPoints;
    pointRemover.getUnrefimentSet
    (
        facesToRestore,
        localFaces,
        localPoints
    );

    // Undo the changes on the faces that are in error.
    pointRemover.setUnrefinement
    (
        localFaces,
        localPoints,
        meshMod
    );

    // Change the mesh (no inflation)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh.clearOut();
    }

    if (meshRefiner_.overwrite())
    {
        mesh.setInstance(meshRefiner_.oldInstance());
    }

    pointRemover.updateMesh(map);
    meshRefiner_.updateMesh(map, labelList(0));

    return map;
}


// Collect all faces that are both in candidateFaces and in the set.
// If coupled face also collects the coupled face.
Foam::labelList Foam::autoLayerDriver::collectFaces
(
    const labelList& candidateFaces,
    const labelHashSet& set
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    // Has face been selected?
    boolList selected(mesh.nFaces(), false);

    forAll(candidateFaces, i)
    {
        label faceI = candidateFaces[i];

        if (set.found(faceI))
        {
            selected[faceI] = true;
        }
    }
    syncTools::syncFaceList
    (
        mesh,
        selected,
        orEqOp<bool>(),     // combine operator
        false               // separation
    );

    labelList selectedFaces(findIndices(selected, true));

    return selectedFaces;
}


// Pick up faces of cells of faces in set.
Foam::labelList Foam::autoLayerDriver::growFaceCellFace
(
    const labelHashSet& set
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    boolList selected(mesh.nFaces(), false);

    forAllConstIter(faceSet, set, iter)
    {
        label faceI = iter.key();

        label own = mesh.faceOwner()[faceI];

        const cell& ownFaces = mesh.cells()[own];
        forAll(ownFaces, ownFaceI)
        {
            selected[ownFaces[ownFaceI]] = true;
        }

        if (mesh.isInternalFace(faceI))
        {
            label nbr = mesh.faceNeighbour()[faceI];

            const cell& nbrFaces = mesh.cells()[nbr];
            forAll(nbrFaces, nbrFaceI)
            {
                selected[nbrFaces[nbrFaceI]] = true;
            }
        }
    }
    syncTools::syncFaceList
    (
        mesh,
        selected,
        orEqOp<bool>(),     // combine operator
        false               // separation
    );
    return findIndices(selected, true);
}


// Remove points not used by any face or points used by only two faces where
// the edges are in line
Foam::label Foam::autoLayerDriver::mergeEdgesUndo
(
    const scalar minCos,
    const dictionary& motionDict
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl
        << "Merging all points on surface that" << nl
        << "- are used by only two boundary faces and" << nl
        << "- make an angle with a cosine of more than " << minCos
        << "." << nl << endl;

    // Point removal analysis engine with undo
    removePoints pointRemover(mesh, true);

    // Count usage of points
    boolList pointCanBeDeleted;
    label nRemove = pointRemover.countPointUsage(minCos, pointCanBeDeleted);

    if (nRemove > 0)
    {
        Info<< "Removing " << nRemove
            << " straight edge points ..." << nl << endl;

        // Remove points
        // ~~~~~~~~~~~~~

        doRemovePoints(pointRemover, pointCanBeDeleted);


        for (label iteration = 0; iteration < 100; iteration++)
        {
            Info<< nl
                << "Undo iteration " << iteration << nl
                << "----------------" << endl;


            // Check mesh for errors
            // ~~~~~~~~~~~~~~~~~~~~~

            faceSet errorFaces
            (
                mesh,
                "errorFaces",
                mesh.nFaces()-mesh.nInternalFaces()
            );
            bool hasErrors = motionSmoother::checkMesh
            (
                false,  // report
                mesh,
                motionDict,
                errorFaces
            );
            //if (checkEdgeConnectivity)
            //{
            //    Info<< "Checking edge-face connectivity (duplicate faces"
            //        << " or non-consecutive shared vertices)" << endl;
            //
            //    label nOldSize = errorFaces.size();
            //
            //    hasErrors =
            //        mesh.checkFaceFaces
            //        (
            //            false,
            //            &errorFaces
            //        )
            //     || hasErrors;
            //
            //    Info<< "Detected additional "
            //        << returnReduce(errorFaces.size()-nOldSize,sumOp<label>())
            //        << " faces with illegal face-face connectivity"
            //        << endl;
            //}

            if (!hasErrors)
            {
                break;
            }

            if (debug)
            {
                Pout<< "**Writing all faces in error to faceSet "
                    << errorFaces.objectPath() << nl << endl;
                errorFaces.write();
            }

            labelList masterErrorFaces
            (
                collectFaces
                (
                    pointRemover.savedFaceLabels(),
                    errorFaces
                )
            );

            label n = returnReduce(masterErrorFaces.size(), sumOp<label>());

            Info<< "Detected " << n
                << " error faces on boundaries that have been merged."
                << " These will be restored to their original faces." << nl
                << endl;

            if (n == 0)
            {
                if (hasErrors)
                {
                    Info<< "Detected "
                        << returnReduce(errorFaces.size(), sumOp<label>())
                        << " error faces in mesh."
                        << " Restoring neighbours of faces in error." << nl
                        << endl;

                    labelList expandedErrorFaces
                    (
                        growFaceCellFace
                        (
                            errorFaces
                        )
                    );

                    doRestorePoints(pointRemover, expandedErrorFaces);
                }

                break;
            }

            doRestorePoints(pointRemover, masterErrorFaces);
        }

        if (debug)
        {
            Pout<< "Writing merged-edges mesh to time "
                << meshRefiner_.timeName() << nl << endl;
            mesh.write();
        }
    }
    else
    {
        Info<< "No straight edges simplified and no points removed ..." << endl;
    }

    return nRemove;
}


// For debugging: Dump displacement to .obj files
void Foam::autoLayerDriver::dumpDisplacement
(
    const fileName& prefix,
    const indirectPrimitivePatch& pp,
    const vectorField& patchDisp,
    const List<extrudeMode>& extrudeStatus
)
{
    OFstream dispStr(prefix + "_disp.obj");
    Info<< "Writing all displacements to " << dispStr.name() << nl << endl;

    label vertI = 0;

    forAll(patchDisp, patchPointI)
    {
        const point& pt = pp.localPoints()[patchPointI];

        meshTools::writeOBJ(dispStr, pt); vertI++;
        meshTools::writeOBJ(dispStr, pt + patchDisp[patchPointI]); vertI++;

        dispStr << "l " << vertI-1 << ' ' << vertI << nl;
    }


    OFstream illStr(prefix + "_illegal.obj");
    Info<< "Writing invalid displacements to " << illStr.name() << nl << endl;

    vertI = 0;

    forAll(patchDisp, patchPointI)
    {
        if (extrudeStatus[patchPointI] != EXTRUDE)
        {
            const point& pt = pp.localPoints()[patchPointI];

            meshTools::writeOBJ(illStr, pt); vertI++;
            meshTools::writeOBJ(illStr, pt + patchDisp[patchPointI]); vertI++;

            illStr << "l " << vertI-1 << ' ' << vertI << nl;
        }
    }
}


// Check that primitivePatch is not multiply connected. Collect non-manifold
// points in pointSet.
void Foam::autoLayerDriver::checkManifold
(
    const indirectPrimitivePatch& fp,
    pointSet& nonManifoldPoints
)
{
    // Check for non-manifold points (surface pinched at point)
    fp.checkPointManifold(false, &nonManifoldPoints);

    // Check for edge-faces (surface pinched at edge)
    const labelListList& edgeFaces = fp.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() > 2)
        {
            const edge& e = fp.edges()[edgeI];

            nonManifoldPoints.insert(fp.meshPoints()[e[0]]);
            nonManifoldPoints.insert(fp.meshPoints()[e[1]]);
        }
    }
}


void Foam::autoLayerDriver::checkMeshManifold() const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Checking mesh manifoldness ..." << endl;

    // Get all outside faces
    labelList outsideFaces(mesh.nFaces() - mesh.nInternalFaces());

    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
         outsideFaces[faceI - mesh.nInternalFaces()] = faceI;
    }

    pointSet nonManifoldPoints
    (
        mesh,
        "nonManifoldPoints",
        mesh.nPoints() / 100
    );

    // Build primitivePatch out of faces and check it for problems.
    checkManifold
    (
        indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), outsideFaces),
            mesh.points()
        ),
        nonManifoldPoints
    );

    label nNonManif = returnReduce(nonManifoldPoints.size(), sumOp<label>());

    if (nNonManif > 0)
    {
        Info<< "Outside of mesh is multiply connected across edges or"
            << " points." << nl
            << "This is not a fatal error but might cause some unexpected"
            << " behaviour." << nl
            << "Writing " << nNonManif
            << " points where this happens to pointSet "
            << nonManifoldPoints.name() << endl;

        nonManifoldPoints.write();
    }
    Info<< endl;
}



// Unset extrusion on point. Returns true if anything unset.
bool Foam::autoLayerDriver::unmarkExtrusion
(
    const label patchPointI,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
)
{
    if (extrudeStatus[patchPointI] == EXTRUDE)
    {
        extrudeStatus[patchPointI] = NOEXTRUDE;
        patchNLayers[patchPointI] = 0;
        patchDisp[patchPointI] = vector::zero;
        return true;
    }
    else if (extrudeStatus[patchPointI] == EXTRUDEREMOVE)
    {
        extrudeStatus[patchPointI] = NOEXTRUDE;
        patchNLayers[patchPointI] = 0;
        patchDisp[patchPointI] = vector::zero;
        return true;
    }
    else
    {
        return false;
    }
}


// Unset extrusion on face. Returns true if anything unset.
bool Foam::autoLayerDriver::unmarkExtrusion
(
    const face& localFace,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
)
{
    bool unextruded = false;

    forAll(localFace, fp)
    {
        if
        (
            unmarkExtrusion
            (
                localFace[fp],
                patchDisp,
                patchNLayers,
                extrudeStatus
            )
        )
        {
            unextruded = true;
        }
    }
    return unextruded;
}


// No extrusion at non-manifold points.
void Foam::autoLayerDriver::handleNonManifolds
(
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Handling non-manifold points ..." << endl;

    // Detect non-manifold points
    Info<< nl << "Checking patch manifoldness ..." << endl;

    pointSet nonManifoldPoints(mesh, "nonManifoldPoints", pp.nPoints());

    // 1. Local check
    checkManifold(pp, nonManifoldPoints);

    label nNonManif = returnReduce(nonManifoldPoints.size(), sumOp<label>());

    Info<< "Outside of local patch is multiply connected across edges or"
        << " points at " << nNonManif << " points." << endl;


    const labelList& meshPoints = pp.meshPoints();


    // 2. Parallel check
    // For now disable extrusion at any shared edges.
    const labelHashSet sharedEdgeSet(mesh.globalData().sharedEdgeLabels());

    forAll(pp.edges(), edgeI)
    {
        if (sharedEdgeSet.found(meshEdges[edgeI]))
        {
            const edge& e = mesh.edges()[meshEdges[edgeI]];

            Pout<< "Disabling extrusion at edge "
                << mesh.points()[e[0]] << mesh.points()[e[1]]
                << " since it is non-manifold across coupled patches."
                << endl;

            nonManifoldPoints.insert(e[0]);
            nonManifoldPoints.insert(e[1]);
        }
    }

    // 3b. extrusion can produce multiple faces between 2 cells
    // across processor boundary
    // This occurs when a coupled face shares more than 1 edge with a
    // non-coupled boundary face.
    // This is now correctly handled by addPatchCellLayer in that it
    // extrudes a single face from the stringed up edges.


    nNonManif = returnReduce(nonManifoldPoints.size(), sumOp<label>());

    if (nNonManif > 0)
    {
        Info<< "Outside of patches is multiply connected across edges or"
            << " points at " << nNonManif << " points." << nl
            << "Writing " << nNonManif
            << " points where this happens to pointSet "
            << nonManifoldPoints.name() << nl
            << "and setting layer thickness to zero on these points."
            << endl;

        nonManifoldPoints.write();

        forAll(meshPoints, patchPointI)
        {
            if (nonManifoldPoints.found(meshPoints[patchPointI]))
            {
                unmarkExtrusion
                (
                    patchPointI,
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                );
            }
        }
    }

    Info<< "Set displacement to zero for all " << nNonManif
        << " non-manifold points" << endl;
}


// Parallel feature edge detection. Assumes non-manifold edges already handled.
void Foam::autoLayerDriver::handleFeatureAngle
(
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const scalar minCos,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Handling feature edges ..." << endl;

    if (minCos < 1-SMALL)
    {
        // Normal component of normals of connected faces.
        vectorField edgeNormal(mesh.nEdges(), wallPoint::greatPoint);

        const labelListList& edgeFaces = pp.edgeFaces();

        forAll(edgeFaces, edgeI)
        {
            const labelList& eFaces = pp.edgeFaces()[edgeI];

            label meshEdgeI = meshEdges[edgeI];

            forAll(eFaces, i)
            {
                nomalsCombine()
                (
                    edgeNormal[meshEdgeI],
                    pp.faceNormals()[eFaces[i]]
                );
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            edgeNormal,
            nomalsCombine(),
            wallPoint::greatPoint,  // null value
            false                   // no separation
        );

        label vertI = 0;
        autoPtr<OFstream> str;
        if (debug)
        {
            str.reset(new OFstream(mesh.time().path()/"featureEdges.obj"));
            Info<< "Writing feature edges to " << str().name() << endl;
        }

        label nFeats = 0;

        // Now on coupled edges the edgeNormal will have been truncated and
        // only be still be the old value where two faces have the same normal
        forAll(edgeFaces, edgeI)
        {
            const labelList& eFaces = pp.edgeFaces()[edgeI];

            label meshEdgeI = meshEdges[edgeI];

            const vector& n = edgeNormal[meshEdgeI];

            if (n != wallPoint::greatPoint)
            {
                scalar cos = n & pp.faceNormals()[eFaces[0]];

                if (cos < minCos)
                {
                    const edge& e = pp.edges()[edgeI];

                    unmarkExtrusion
                    (
                        e[0],
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    );
                    unmarkExtrusion
                    (
                        e[1],
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    );

                    nFeats++;

                    if (str.valid())
                    {
                        meshTools::writeOBJ(str(), pp.localPoints()[e[0]]);
                        vertI++;
                        meshTools::writeOBJ(str(), pp.localPoints()[e[1]]);
                        vertI++;
                        str()<< "l " << vertI-1 << ' ' << vertI << nl;
                    }
                }
            }
        }

        Info<< "Set displacement to zero for points on "
            << returnReduce(nFeats, sumOp<label>())
            << " feature edges" << endl;
    }
}


// No extrusion on cells with warped faces. Calculates the thickness of the
// layer and compares it to the space the warped face takes up. Disables
// extrusion if layer thickness is more than faceRatio of the thickness of
// the face.
void Foam::autoLayerDriver::handleWarpedFaces
(
    const indirectPrimitivePatch& pp,
    const scalar faceRatio,
    const scalar edge0Len,
    const labelList& cellLevel,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Handling cells with warped patch faces ..." << nl;

    const pointField& points = mesh.points();

    label nWarpedFaces = 0;

    forAll(pp, i)
    {
        const face& f = pp[i];

        if (f.size() > 3)
        {
            label faceI = pp.addressing()[i];

            label ownLevel = cellLevel[mesh.faceOwner()[faceI]];
            scalar edgeLen = edge0Len/(1<<ownLevel);

            // Normal distance to face centre plane
            const point& fc = mesh.faceCentres()[faceI];
            const vector& fn = pp.faceNormals()[i];

            scalarField vProj(f.size());

            forAll(f, fp)
            {
                vector n = points[f[fp]] - fc;
                vProj[fp] = (n & fn);
            }

            // Get normal 'span' of face
            scalar minVal = min(vProj);
            scalar maxVal = max(vProj);

            if ((maxVal - minVal) > faceRatio * edgeLen)
            {
                if
                (
                    unmarkExtrusion
                    (
                        pp.localFaces()[i],
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nWarpedFaces++;
                }
            }
        }
    }

    Info<< "Set displacement to zero on "
        << returnReduce(nWarpedFaces, sumOp<label>())
        << " warped faces since layer would be > " << faceRatio
        << " of the size of the bounding box." << endl;
}



//// No extrusion on cells with multiple patch faces. There ususally is a reason
//// why combinePatchFaces hasn't succeeded.
//void Foam::autoLayerDriver::handleMultiplePatchFaces
//(
//    const indirectPrimitivePatch& pp,
//    pointField& patchDisp,
//    labelList& patchNLayers,
//    List<extrudeMode>& extrudeStatus
//) const
//{
//    const fvMesh& mesh = meshRefiner_.mesh();
//
//    Info<< nl << "Handling cells with multiple patch faces ..." << nl;
//
//    const labelListList& pointFaces = pp.pointFaces();
//
//    // Cells that should not get an extrusion layer
//    cellSet multiPatchCells(mesh, "multiPatchCells", pp.size());
//
//    // Detect points that use multiple faces on same cell.
//    forAll(pointFaces, patchPointI)
//    {
//        const labelList& pFaces = pointFaces[patchPointI];
//
//        labelHashSet pointCells(pFaces.size());
//
//        forAll(pFaces, i)
//        {
//            label cellI = mesh.faceOwner()[pp.addressing()[pFaces[i]]];
//
//            if (!pointCells.insert(cellI))
//            {
//                // Second or more occurrence of cell so cell has two or more
//                // pp faces connected to this point.
//                multiPatchCells.insert(cellI);
//            }
//        }
//    }
//
//    label nMultiPatchCells = returnReduce
//    (
//        multiPatchCells.size(),
//        sumOp<label>()
//    );
//
//    Info<< "Detected " << nMultiPatchCells
//        << " cells with multiple (connected) patch faces." << endl;
//
//    label nChanged = 0;
//
//    if (nMultiPatchCells > 0)
//    {
//        Info<< "Writing " << nMultiPatchCells
//            << " cells with multiple (connected) patch faces to cellSet "
//            << multiPatchCells.objectPath() << endl;
//        multiPatchCells.write();
//
//
//        // Go through all points and remove extrusion on any cell in
//        // multiPatchCells
//        // (has to be done in separate loop since having one point on
//        // multipatches has to reset extrusion on all points of cell)
//
//        forAll(pointFaces, patchPointI)
//        {
//            if (extrudeStatus[patchPointI] != NOEXTRUDE)
//            {
//                const labelList& pFaces = pointFaces[patchPointI];
//
//                forAll(pFaces, i)
//                {
//                    label cellI =
//                        mesh.faceOwner()[pp.addressing()[pFaces[i]]];
//
//                    if (multiPatchCells.found(cellI))
//                    {
//                        if
//                        (
//                            unmarkExtrusion
//                            (
//                                patchPointI,
//                                patchDisp,
//                                patchNLayers,
//                                extrudeStatus
//                            )
//                        )
//                        {
//                            nChanged++;
//                        }
//                    }
//                }
//            }
//        }
//
//        reduce(nChanged, sumOp<label>());
//    }
//
//    Info<< "Prevented extrusion on " << nChanged
//        << " points due to multiple patch faces." << nl << endl;
//}


// No extrusion on faces with differing number of layers for points
void Foam::autoLayerDriver::setNumLayers
(
    const labelList& patchToNLayers,
    const labelList& patchIDs,
    const indirectPrimitivePatch& pp,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Handling points with inconsistent layer specification ..."
        << endl;

    // Get for every point (really only nessecary on patch external points)
    // the max and min of any patch faces using it.
    labelList maxLayers(patchNLayers.size(), labelMin);
    labelList minLayers(patchNLayers.size(), labelMax);

    forAll(patchIDs, i)
    {
        label patchI = patchIDs[i];

        const labelList& meshPoints = mesh.boundaryMesh()[patchI].meshPoints();

        label wantedLayers = patchToNLayers[patchI];

        forAll(meshPoints, patchPointI)
        {
            label ppPointI = pp.meshPointMap()[meshPoints[patchPointI]];

            maxLayers[ppPointI] = max(wantedLayers, maxLayers[ppPointI]);
            minLayers[ppPointI] = min(wantedLayers, minLayers[ppPointI]);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        maxLayers,
        maxEqOp<label>(),
        labelMin,           // null value
        false               // no separation
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        minLayers,
        minEqOp<label>(),
        labelMax,           // null value
        false               // no separation
    );

    // Unmark any point with different min and max
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //label nConflicts = 0;

    forAll(maxLayers, i)
    {
        if (maxLayers[i] == labelMin || minLayers[i] == labelMax)
        {
            FatalErrorIn("setNumLayers(..)")
                << "Patchpoint:" << i << " coord:" << pp.localPoints()[i]
                << " maxLayers:" << maxLayers
                << " minLayers:" << minLayers
                << abort(FatalError);
        }
        else if (maxLayers[i] == minLayers[i])
        {
            // Ok setting.
            patchNLayers[i] = maxLayers[i];
        }
        else
        {
            // Inconsistent num layers between patch faces using point
            //if
            //(
            //    unmarkExtrusion
            //    (
            //        i,
            //        patchDisp,
            //        patchNLayers,
            //        extrudeStatus
            //    )
            //)
            //{
            //    nConflicts++;
            //}
            patchNLayers[i] = maxLayers[i];
        }
    }

    //reduce(nConflicts, sumOp<label>());
    //
    //Info<< "Set displacement to zero for " << nConflicts
    //    << " points due to points being on multiple regions"
    //    << " with inconsistent nLayers specification." << endl;
}


// Grow no-extrusion layer
void Foam::autoLayerDriver::growNoExtrusion
(
    const indirectPrimitivePatch& pp,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
)
{
    Info<< nl << "Growing non-extrusion points by one layer ..." << endl;

    List<extrudeMode> grownExtrudeStatus(extrudeStatus);

    const faceList& localFaces = pp.localFaces();

    label nGrown = 0;

    forAll(localFaces, faceI)
    {
        const face& f = localFaces[faceI];

        bool hasSqueeze = false;
        forAll(f, fp)
        {
            if (extrudeStatus[f[fp]] == NOEXTRUDE)
            {
                hasSqueeze = true;
                break;
            }
        }

        if (hasSqueeze)
        {
            // Squeeze all points of face
            forAll(f, fp)
            {
                if
                (
                    extrudeStatus[f[fp]] == NOEXTRUDE
                 && grownExtrudeStatus[f[fp]] != NOEXTRUDE
                )
                {
                    grownExtrudeStatus[f[fp]] = NOEXTRUDE;
                    nGrown++;
                }
            }
        }
    }

    extrudeStatus.transfer(grownExtrudeStatus);

    forAll(extrudeStatus, patchPointI)
    {
        if (extrudeStatus[patchPointI] == NOEXTRUDE)
        {
            patchDisp[patchPointI] = vector::zero;
            patchNLayers[patchPointI] = 0;
        }
    }

    reduce(nGrown, sumOp<label>());

    Info<< "Set displacement to zero for an additional " << nGrown
        << " points." << endl;
}


void Foam::autoLayerDriver::calculateLayerThickness
(
    const indirectPrimitivePatch& pp,
    const labelList& patchIDs,
    const scalarField& patchExpansionRatio,

    const bool relativeSizes,
    const scalarField& patchFinalLayerThickness,
    const scalarField& patchMinThickness,

    const labelList& cellLevel,
    const labelList& patchNLayers,
    const scalar edge0Len,

    scalarField& thickness,
    scalarField& minThickness,
    scalarField& expansionRatio
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();


    // Rework patch-wise layer parameters into minimum per point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Reuse input fields
    expansionRatio.setSize(pp.nPoints());
    expansionRatio = GREAT;
    thickness.setSize(pp.nPoints());
    thickness = GREAT;
    minThickness.setSize(pp.nPoints());
    minThickness = GREAT;

    forAll(patchIDs, i)
    {
        label patchI = patchIDs[i];

        const labelList& meshPoints = patches[patchI].meshPoints();

        forAll(meshPoints, patchPointI)
        {
            label ppPointI = pp.meshPointMap()[meshPoints[patchPointI]];

            expansionRatio[ppPointI] = min
            (
                expansionRatio[ppPointI],
                patchExpansionRatio[patchI]
            );
            thickness[ppPointI] = min
            (
                thickness[ppPointI],
                patchFinalLayerThickness[patchI]
            );
            minThickness[ppPointI] = min
            (
                minThickness[ppPointI],
                patchMinThickness[patchI]
            );
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        expansionRatio,
        minEqOp<scalar>(),
        GREAT,              // null value
        false               // no separation
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        thickness,
        minEqOp<scalar>(),
        GREAT,              // null value
        false               // no separation
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        minThickness,
        minEqOp<scalar>(),
        GREAT,              // null value
        false               // no separation
    );


    // Now the thicknesses are set according to the minimum of connected
    // patches.


    // Rework relative thickness into absolute
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // by multiplying with the internal cell size.

    if (relativeSizes)
    {
        if (min(patchMinThickness) < 0 || max(patchMinThickness) > 2)
        {
            FatalErrorIn("calculateLayerThickness(..)")
                << "Thickness should be factor of local undistorted cell size."
                << " Valid values are [0..2]." << nl
                << " minThickness:" << patchMinThickness
                << exit(FatalError);
        }


        // Determine per point the max cell level of connected cells
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        labelList maxPointLevel(pp.nPoints(), labelMin);

        forAll(pp, i)
        {
            label ownLevel = cellLevel[mesh.faceOwner()[pp.addressing()[i]]];

            const face& f = pp.localFaces()[i];

            forAll(f, fp)
            {
                maxPointLevel[f[fp]] = max(maxPointLevel[f[fp]], ownLevel);
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            maxPointLevel,
            maxEqOp<label>(),
            labelMin,           // null value
            false               // no separation
        );


        forAll(maxPointLevel, pointI)
        {
            // Find undistorted edge size for this level.
            scalar edgeLen = edge0Len/(1<<maxPointLevel[pointI]);
            thickness[pointI] *= edgeLen;
            minThickness[pointI] *= edgeLen;
        }
    }



    // Rework thickness (of final layer) into overall thickness of all layers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(thickness, pointI)
    {
        // Calculate layer thickness based on expansion ratio
        // and final layer height
        if (expansionRatio[pointI] == 1)
        {
            thickness[pointI] *= patchNLayers[pointI];
        }
        else
        {

            scalar invExpansion = 1.0 / expansionRatio[pointI];
            label nLay = patchNLayers[pointI];
            thickness[pointI] *=
                (1.0 - pow(invExpansion, nLay))
              / (1.0 - invExpansion);
        }
    }


    //Info<< "calculateLayerThickness : min:" << gMin(thickness)
    //    << " max:" << gMax(thickness) << endl;
}


// Synchronize displacement among coupled patches.
void Foam::autoLayerDriver::syncPatchDisplacement
(
    const motionSmoother& meshMover,
    const scalarField& minThickness,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const labelList& meshPoints = meshMover.patch().meshPoints();

    label nChangedTotal = 0;

    while (true)
    {
        label nChanged = 0;

        // Sync displacement (by taking min)
        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            patchDisp,
            minEqOp<vector>(),
            wallPoint::greatPoint,          // null value
            false                           // no separation
        );

        // Unmark if displacement too small
        forAll(patchDisp, i)
        {
            if (mag(patchDisp[i]) < minThickness[i])
            {
                if
                (
                    unmarkExtrusion
                    (
                        i,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nChanged++;
                }
            }
        }

        labelList syncPatchNLayers(patchNLayers);

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            syncPatchNLayers,
            minEqOp<label>(),
            labelMax,           // null value
            false               // no separation
        );

        // Reset if differs
        forAll(syncPatchNLayers, i)
        {
            if (syncPatchNLayers[i] != patchNLayers[i])
            {
                if
                (
                    unmarkExtrusion
                    (
                        i,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nChanged++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            syncPatchNLayers,
            maxEqOp<label>(),
            labelMin,           // null value
            false               // no separation
        );

        // Reset if differs
        forAll(syncPatchNLayers, i)
        {
            if (syncPatchNLayers[i] != patchNLayers[i])
            {
                if
                (
                    unmarkExtrusion
                    (
                        i,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nChanged++;
                }
            }
        }
        nChangedTotal += nChanged;

        if (!returnReduce(nChanged, sumOp<label>()))
        {
            break;
        }
    }

    Info<< "Prevented extrusion on "
        << returnReduce(nChangedTotal, sumOp<label>())
        << " coupled patch points during syncPatchDisplacement." << endl;
}


// Calculate displacement vector for all patch points. Uses pointNormal.
// Checks that displaced patch point would be visible from all centres
// of the faces using it.
// extrudeStatus is both input and output and gives the status of each
// patch point.
void Foam::autoLayerDriver::getPatchDisplacement
(
    const motionSmoother& meshMover,
    const scalarField& thickness,
    const scalarField& minThickness,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    Info<< nl << "Determining displacement for added points"
        << " according to pointNormal ..." << endl;

    const fvMesh& mesh = meshRefiner_.mesh();
    const indirectPrimitivePatch& pp = meshMover.patch();
    const vectorField& faceNormals = pp.faceNormals();
    const labelListList& pointFaces = pp.pointFaces();
    const pointField& localPoints = pp.localPoints();
    const labelList& meshPoints = pp.meshPoints();

    // Determine pointNormal
    // ~~~~~~~~~~~~~~~~~~~~~

    pointField pointNormals(pp.nPoints(), vector::zero);
    {
        labelList nPointFaces(pp.nPoints(), 0);

        forAll(faceNormals, faceI)
        {
            const face& f = pp.localFaces()[faceI];

            forAll(f, fp)
            {
                pointNormals[f[fp]] += faceNormals[faceI];
                nPointFaces[f[fp]] ++;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            pointNormals,
            plusEqOp<vector>(),
            vector::zero,       // null value
            false               // no separation
        );

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            nPointFaces,
            plusEqOp<label>(),
            0,                  // null value
            false               // no separation
        );

        forAll(pointNormals, i)
        {
            pointNormals[i] /= nPointFaces[i];
        }
    }


    // Determine local length scale on patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Start off from same thickness everywhere (except where no extrusion)
    patchDisp = thickness*pointNormals;

    // Check if no extrude possible.
    forAll(pointNormals, patchPointI)
    {
        label meshPointI = pp.meshPoints()[patchPointI];

        if (extrudeStatus[patchPointI] == NOEXTRUDE)
        {
            // Do not use unmarkExtrusion; forcibly set to zero extrusion.
            patchNLayers[patchPointI] = 0;
            patchDisp[patchPointI] = vector::zero;
        }
        else
        {
            // Get normal
            const vector& n = pointNormals[patchPointI];

            if (!meshTools::visNormal(n, faceNormals, pointFaces[patchPointI]))
            {
                Pout<< "No valid normal for point " << meshPointI
                    << ' ' << pp.points()[meshPointI]
                    << "; setting displacement to " << patchDisp[patchPointI]
                    << endl;

                extrudeStatus[patchPointI] = EXTRUDEREMOVE;
            }
        }
    }

    // At illegal points make displacement average of new neighbour positions
    forAll(extrudeStatus, patchPointI)
    {
        if (extrudeStatus[patchPointI] == EXTRUDEREMOVE)
        {
            point avg(vector::zero);
            label nPoints = 0;

            const labelList& pEdges = pp.pointEdges()[patchPointI];

            forAll(pEdges, i)
            {
                label edgeI = pEdges[i];

                label otherPointI = pp.edges()[edgeI].otherVertex(patchPointI);

                if (extrudeStatus[otherPointI] != NOEXTRUDE)
                {
                    avg += localPoints[otherPointI] + patchDisp[otherPointI];
                    nPoints++;
                }
            }

            if (nPoints > 0)
            {
                Pout<< "Displacement at illegal point "
                    << localPoints[patchPointI]
                    << " set to " << (avg / nPoints - localPoints[patchPointI])
                    << endl;

                patchDisp[patchPointI] =
                    avg / nPoints
                  - localPoints[patchPointI];
            }
        }
    }

    // Make sure displacement is equal on both sides of coupled patches.
    syncPatchDisplacement
    (
        meshMover,
        minThickness,
        patchDisp,
        patchNLayers,
        extrudeStatus
    );

    Info<< endl;
}


// Truncates displacement
// - for all patchFaces in the faceset displacement gets set to zero
// - all displacement < minThickness gets set to zero
Foam::label Foam::autoLayerDriver::truncateDisplacement
(
    const motionSmoother& meshMover,
    const scalarField& minThickness,
    const faceSet& illegalPatchFaces,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const polyMesh& mesh = meshMover.mesh();
    const indirectPrimitivePatch& pp = meshMover.patch();

    label nChanged = 0;

    const Map<label>& meshPointMap = pp.meshPointMap();

    forAllConstIter(faceSet, illegalPatchFaces, iter)
    {
        label faceI = iter.key();

        if (mesh.isInternalFace(faceI))
        {
            FatalErrorIn("truncateDisplacement(..)")
                << "Faceset " << illegalPatchFaces.name()
                << " contains internal face " << faceI << nl
                << "It should only contain patch faces" << abort(FatalError);
        }

        const face& f = mesh.faces()[faceI];


        forAll(f, fp)
        {
            if (meshPointMap.found(f[fp]))
            {
                label patchPointI = meshPointMap[f[fp]];

                if (extrudeStatus[patchPointI] != NOEXTRUDE)
                {
                    unmarkExtrusion
                    (
                        patchPointI,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    );
                    nChanged++;
                }
            }
        }
    }

    forAll(patchDisp, patchPointI)
    {
        if (mag(patchDisp[patchPointI]) < minThickness[patchPointI])
        {
            if
            (
                unmarkExtrusion
                (
                    patchPointI,
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                )
            )
            {
                nChanged++;
            }
        }
        else if (extrudeStatus[patchPointI] == NOEXTRUDE)
        {
            // Make sure displacement is 0. Should already be so but ...
            patchDisp[patchPointI] = vector::zero;
            patchNLayers[patchPointI] = 0;
        }
    }


    const faceList& localFaces = pp.localFaces();

    while (true)
    {
        syncPatchDisplacement
        (
            meshMover,
            minThickness,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Make sure that a face doesn't have two non-consecutive areas
        // not extruded (e.g. quad where vertex 0 and 2 are not extruded
        // but 1 and 3 are) since this gives topological errors.

        label nPinched = 0;

        forAll(localFaces, i)
        {
            const face& localF = localFaces[i];

            // Count number of transitions from unsnapped to snapped.
            label nTrans = 0;

            extrudeMode prevMode = extrudeStatus[localF.prevLabel(0)];

            forAll(localF, fp)
            {
                extrudeMode fpMode = extrudeStatus[localF[fp]];

                if (prevMode == NOEXTRUDE && fpMode != NOEXTRUDE)
                {
                    nTrans++;
                }
                prevMode = fpMode;
            }

            if (nTrans > 1)
            {
                // Multiple pinches. Reset whole face as unextruded.
                if
                (
                    unmarkExtrusion
                    (
                        localF,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nPinched++;
                    nChanged++;
                }
            }
        }

        reduce(nPinched, sumOp<label>());

        Info<< "truncateDisplacement : Unextruded " << nPinched
            << " faces due to non-consecutive vertices being extruded." << endl;


        // Make sure that a face has consistent number of layers for all
        // its vertices.

        label nDiffering = 0;

        //forAll(localFaces, i)
        //{
        //    const face& localF = localFaces[i];
        //
        //    label numLayers = -1;
        //
        //    forAll(localF, fp)
        //    {
        //        if (patchNLayers[localF[fp]] > 0)
        //        {
        //            if (numLayers == -1)
        //            {
        //                numLayers = patchNLayers[localF[fp]];
        //            }
        //            else if (numLayers != patchNLayers[localF[fp]])
        //            {
        //                // Differing number of layers
        //                if
        //                (
        //                    unmarkExtrusion
        //                    (
        //                        localF,
        //                        patchDisp,
        //                        patchNLayers,
        //                        extrudeStatus
        //                    )
        //                )
        //                {
        //                    nDiffering++;
        //                    nChanged++;
        //                }
        //                break;
        //            }
        //        }
        //    }
        //}
        //
        //reduce(nDiffering, sumOp<label>());
        //
        //Info<< "truncateDisplacement : Unextruded " << nDiffering
        //    << " faces due to having differing number of layers." << endl;

        if (nPinched+nDiffering == 0)
        {
            break;
        }
    }

    return nChanged;
}


// Setup layer information (at points and faces) to modify mesh topology in
// regions where layer mesh terminates.
void Foam::autoLayerDriver::setupLayerInfoTruncation
(
    const motionSmoother& meshMover,
    const labelList& patchNLayers,
    const List<extrudeMode>& extrudeStatus,
    const label nBufferCellsNoExtrude,
    labelList& nPatchPointLayers,
    labelList& nPatchFaceLayers
) const
{
    Info<< nl << "Setting up information for layer truncation ..." << endl;

    const indirectPrimitivePatch& pp = meshMover.patch();
    const polyMesh& mesh = meshMover.mesh();

    if (nBufferCellsNoExtrude < 0)
    {
        Info<< nl << "Performing no layer truncation."
            << " nBufferCellsNoExtrude set to less than 0  ..." << endl;

        // Face layers if any point get extruded
        forAll(pp.localFaces(), patchFaceI)
        {
            const face& f = pp.localFaces()[patchFaceI];

            forAll(f, fp)
            {
                if (patchNLayers[f[fp]] > 0)
                {
                    nPatchFaceLayers[patchFaceI] = patchNLayers[f[fp]];
                    break;
                }
            }
        }
        nPatchPointLayers = patchNLayers;
    }
    else
    {
        // Determine max point layers per face.
        labelList maxLevel(pp.size(), 0);

        forAll(pp.localFaces(), patchFaceI)
        {
            const face& f = pp.localFaces()[patchFaceI];

            // find patch faces where layer terminates (i.e contains extrude
            // and noextrude points).

            bool noExtrude = false;
            label mLevel = 0;

            forAll(f, fp)
            {
                if (extrudeStatus[f[fp]] == NOEXTRUDE)
                {
                    noExtrude = true;
                }
                mLevel = max(mLevel, patchNLayers[f[fp]]);
            }

            if (mLevel > 0)
            {
                // So one of the points is extruded. Check if all are extruded
                // or is a mix.

                if (noExtrude)
                {
                    nPatchFaceLayers[patchFaceI] = 1;
                    maxLevel[patchFaceI] = mLevel;
                }
                else
                {
                    maxLevel[patchFaceI] = mLevel;
                }
            }
        }

        // We have the seed faces (faces with nPatchFaceLayers != maxLevel)
        // Now do a meshwave across the patch where we pick up neighbours
        // of seed faces.
        // Note: quite inefficient. Could probably be coded better.

        const labelListList& pointFaces = pp.pointFaces();

        label nLevels = gMax(patchNLayers);

        // flag neighbouring patch faces with number of layers to grow
        for (label ilevel = 1; ilevel < nLevels; ilevel++)
        {
            label nBuffer;

            if (ilevel == 1)
            {
                nBuffer = nBufferCellsNoExtrude - 1;
            }
            else
            {
                nBuffer = nBufferCellsNoExtrude;
            }

            for (label ibuffer = 0; ibuffer < nBuffer + 1; ibuffer++)
            {
                labelList tempCounter(nPatchFaceLayers);

                boolList foundNeighbour(pp.nPoints(), false);

                forAll(pp.meshPoints(), patchPointI)
                {
                    forAll(pointFaces[patchPointI], pointFaceI)
                    {
                        label faceI = pointFaces[patchPointI][pointFaceI];

                        if
                        (
                            nPatchFaceLayers[faceI] != -1
                         && maxLevel[faceI] > 0
                        )
                        {
                            foundNeighbour[patchPointI] = true;
                            break;
                        }
                    }
                }

                syncTools::syncPointList
                (
                    mesh,
                    pp.meshPoints(),
                    foundNeighbour,
                    orEqOp<bool>(),
                    false,              // null value
                    false               // no separation
                );

                forAll(pp.meshPoints(), patchPointI)
                {
                    if (foundNeighbour[patchPointI])
                    {
                        forAll(pointFaces[patchPointI], pointFaceI)
                        {
                            label faceI = pointFaces[patchPointI][pointFaceI];
                            if
                            (
                                nPatchFaceLayers[faceI] == -1
                             && maxLevel[faceI] > 0
                             && ilevel < maxLevel[faceI]
                            )
                            {
                                tempCounter[faceI] = ilevel;
                            }
                        }
                    }
                }
                nPatchFaceLayers = tempCounter;
            }
        }

        forAll(pp.localFaces(), patchFaceI)
        {
            if (nPatchFaceLayers[patchFaceI] == -1)
            {
                nPatchFaceLayers[patchFaceI] = maxLevel[patchFaceI];
            }
        }

        forAll(pp.meshPoints(), patchPointI)
        {
            if (extrudeStatus[patchPointI] != NOEXTRUDE)
            {
                forAll(pointFaces[patchPointI], pointFaceI)
                {
                    label face = pointFaces[patchPointI][pointFaceI];
                    nPatchPointLayers[patchPointI] = max
                    (
                        nPatchPointLayers[patchPointI],
                        nPatchFaceLayers[face]
                    );
                }
            }
            else
            {
                nPatchPointLayers[patchPointI] = 0;
            }
        }
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            nPatchPointLayers,
            maxEqOp<label>(),
            0,                  // null value
            false               // no separation
        );
    }
}


// Does any of the cells use a face from faces?
bool Foam::autoLayerDriver::cellsUseFace
(
    const polyMesh& mesh,
    const labelList& cellLabels,
    const labelHashSet& faces
)
{
    forAll(cellLabels, i)
    {
        const cell& cFaces = mesh.cells()[cellLabels[i]];

        forAll(cFaces, cFaceI)
        {
            if (faces.found(cFaces[cFaceI]))
            {
                return true;
            }
        }
    }
    return false;
}


// Checks the newly added cells and locally unmarks points so they
// will not get extruded next time round. Returns global number of unmarked
// points (0 if all was fine)
Foam::label Foam::autoLayerDriver::checkAndUnmark
(
    const addPatchCellLayer& addLayer,
    const dictionary& meshQualityDict,
    const indirectPrimitivePatch& pp,
    const fvMesh& newMesh,

    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
)
{
    // Check the resulting mesh for errors
    Info<< nl << "Checking mesh with layer ..." << endl;
    faceSet wrongFaces(newMesh, "wrongFaces", newMesh.nFaces()/1000);
    motionSmoother::checkMesh(false, newMesh, meshQualityDict, wrongFaces);
    Info<< "Detected " << returnReduce(wrongFaces.size(), sumOp<label>())
        << " illegal faces"
        << " (concave, zero area or negative cell pyramid volume)"
        << endl;

    // Undo local extrusion if
    // - any of the added cells in error

    label nChanged = 0;

    // Get all cells in the layer.
    labelListList addedCells
    (
        addPatchCellLayer::addedCells
        (
            newMesh,
            addLayer.layerFaces()
        )
    );

    // Check if any of the faces in error uses any face of an added cell
    forAll(addedCells, oldPatchFaceI)
    {
        // Get the cells (in newMesh labels) per old patch face (in mesh
        // labels)
        const labelList& fCells = addedCells[oldPatchFaceI];

        if (cellsUseFace(newMesh, fCells, wrongFaces))
        {
            // Unmark points on old mesh
            if
            (
                unmarkExtrusion
                (
                    pp.localFaces()[oldPatchFaceI],
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                )
            )
            {
                nChanged++;
            }
        }
    }

    return returnReduce(nChanged, sumOp<label>());
}


//- Count global number of extruded faces
Foam::label Foam::autoLayerDriver::countExtrusion
(
    const indirectPrimitivePatch& pp,
    const List<extrudeMode>& extrudeStatus
)
{
    // Count number of extruded patch faces
    label nExtruded = 0;
    {
        const faceList& localFaces = pp.localFaces();

        forAll(localFaces, i)
        {
            const face& localFace = localFaces[i];

            forAll(localFace, fp)
            {
                if (extrudeStatus[localFace[fp]] != NOEXTRUDE)
                {
                    nExtruded++;
                    break;
                }
            }
        }
    }

    return returnReduce(nExtruded, sumOp<label>());
}


// Collect layer faces and layer cells into bools for ease of handling
void Foam::autoLayerDriver::getLayerCellsFaces
(
    const polyMesh& mesh,
    const addPatchCellLayer& addLayer,
    boolList& flaggedCells,
    boolList& flaggedFaces
)
{
    flaggedCells.setSize(mesh.nCells());
    flaggedCells = false;
    flaggedFaces.setSize(mesh.nFaces());
    flaggedFaces = false;

    // Mark all faces in the layer
    const labelListList& layerFaces = addLayer.layerFaces();

    // Mark all cells in the layer.
    labelListList addedCells(addPatchCellLayer::addedCells(mesh, layerFaces));

    forAll(addedCells, oldPatchFaceI)
    {
        const labelList& added = addedCells[oldPatchFaceI];

        forAll(added, i)
        {
            flaggedCells[added[i]] = true;
        }
    }

    forAll(layerFaces, oldPatchFaceI)
    {
        const labelList& layer = layerFaces[oldPatchFaceI];

        if (layer.size())
        {
            for (label i = 1; i < layer.size()-1; i++)
            {
                flaggedFaces[layer[i]] = true;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoLayerDriver::autoLayerDriver(meshRefinement& meshRefiner)
:
    meshRefiner_(meshRefiner)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::autoLayerDriver::mergePatchFacesUndo
(
    const layerParameters& layerParams,
    const dictionary& motionDict
)
{
    scalar minCos = Foam::cos
    (
        layerParams.featureAngle()
      * mathematicalConstant::pi/180.0
    );

    scalar concaveCos = Foam::cos
    (
        layerParams.concaveAngle()
      * mathematicalConstant::pi/180.0
    );

    Info<< nl
        << "Merging all faces of a cell" << nl
        << "---------------------------" << nl
        << "    - which are on the same patch" << nl
        << "    - which make an angle < " << layerParams.featureAngle()
        << " degrees"
        << nl
        << "      (cos:" << minCos << ')' << nl
        << "    - as long as the resulting face doesn't become concave"
        << " by more than "
        << layerParams.concaveAngle() << " degrees" << nl
        << "      (0=straight, 180=fully concave)" << nl
        << endl;

    label nChanged = mergePatchFacesUndo(minCos, concaveCos, motionDict);

    nChanged += mergeEdgesUndo(minCos, motionDict);
}


void Foam::autoLayerDriver::addLayers
(
    const layerParameters& layerParams,
    const dictionary& motionDict,
    const labelList& patchIDs,
    const label nAllowableErrors,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    autoPtr<indirectPrimitivePatch> pp
    (
        meshRefinement::makePatch
        (
            mesh,
            patchIDs
        )
    );

    // Construct iterative mesh mover.
    Info<< "Constructing mesh displacer ..." << endl;

    autoPtr<motionSmoother> meshMover
    (
        new motionSmoother
        (
            mesh,
            pp(),
            patchIDs,
            meshRefinement::makeDisplacementField
            (
                pointMesh::New(mesh),
                patchIDs
            ),
            motionDict
        )
    );


    // Point-wise extrusion data
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Displacement for all pp.localPoints.
    vectorField patchDisp(pp().nPoints(), vector::one);

    // Number of layers for all pp.localPoints. Note: only valid if
    // extrudeStatus = EXTRUDE.
    labelList patchNLayers(pp().nPoints(), 0);

    // Whether to add edge for all pp.localPoints.
    List<extrudeMode> extrudeStatus(pp().nPoints(), EXTRUDE);

    {
        // Get number of layer per point from number of layers per patch
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        setNumLayers
        (
            layerParams.numLayers(),    // per patch the num layers
            meshMover().adaptPatchIDs(),// patches that are being moved
            pp,                         // indirectpatch for all faces moving

            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Precalculate mesh edge labels for patch edges
        labelList meshEdges(pp().meshEdges(mesh.edges(), mesh.pointEdges()));

        // Disable extrusion on non-manifold points
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        handleNonManifolds
        (
            pp,
            meshEdges,

            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Disable extrusion on feature angles
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        handleFeatureAngle
        (
            pp,
            meshEdges,
            layerParams.featureAngle()*mathematicalConstant::pi/180.0,

            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Disable extrusion on warped faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Undistorted edge length
        const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
        const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

        handleWarpedFaces
        (
            pp,
            layerParams.maxFaceThicknessRatio(),
            edge0Len,
            cellLevel,

            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        //// Disable extrusion on cells with multiple patch faces
        //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //
        //handleMultiplePatchFaces
        //(
        //    pp,
        //
        //    patchDisp,
        //    patchNLayers,
        //    extrudeStatus
        //);


        // Grow out region of non-extrusion
        for (label i = 0; i < layerParams.nGrow(); i++)
        {
            growNoExtrusion
            (
                pp,
                patchDisp,
                patchNLayers,
                extrudeStatus
            );
        }
    }


    // Undistorted edge length
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

    // Determine (wanted) point-wise layer thickness and expansion ratio
    scalarField thickness(pp().nPoints());
    scalarField minThickness(pp().nPoints());
    scalarField expansionRatio(pp().nPoints());
    calculateLayerThickness
    (
        pp,
        meshMover().adaptPatchIDs(),
        layerParams.expansionRatio(),

        layerParams.relativeSizes(),        // thickness relative to cellsize?
        layerParams.finalLayerThickness(),  // wanted thicknes
        layerParams.minThickness(),         // minimum thickness

        cellLevel,
        patchNLayers,
        edge0Len,

        thickness,
        minThickness,
        expansionRatio
    );


    // Print a bit
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // Find maximum length of a patch name, for a nicer output
        label maxPatchNameLen = 0;
        forAll(meshMover().adaptPatchIDs(), i)
        {
            label patchI = meshMover().adaptPatchIDs()[i];
            word patchName = patches[patchI].name();
            maxPatchNameLen = max(maxPatchNameLen, label(patchName.size()));
        }

        Info<< nl
            << setf(ios_base::left) << setw(maxPatchNameLen) << "patch"
            << setw(0) << " faces    layers avg thickness[m]" << nl
            << setf(ios_base::left) << setw(maxPatchNameLen) << " "
            << setw(0) << "                 near-wall overall" << nl
            << setf(ios_base::left) << setw(maxPatchNameLen) << "-----"
            << setw(0) << " -----    ------ --------- -------" << endl;

        forAll(meshMover().adaptPatchIDs(), i)
        {
            label patchI = meshMover().adaptPatchIDs()[i];

            const labelList& meshPoints = patches[patchI].meshPoints();

            //scalar maxThickness = -VGREAT;
            //scalar minThickness = VGREAT;
            scalar sumThickness = 0;
            scalar sumNearWallThickness = 0;

            forAll(meshPoints, patchPointI)
            {
                label ppPointI = pp().meshPointMap()[meshPoints[patchPointI]];

                //maxThickness = max(maxThickness, thickness[ppPointI]);
                //minThickness = min(minThickness, thickness[ppPointI]);
                sumThickness += thickness[ppPointI];

                label nLay = patchNLayers[ppPointI];
                if (nLay > 0)
                {
                    if (expansionRatio[ppPointI] == 1)
                    {
                        sumNearWallThickness += thickness[ppPointI]/nLay;
                    }
                    else
                    {
                        scalar s =
                            (1.0-pow(expansionRatio[ppPointI], nLay))
                          / (1.0-expansionRatio[ppPointI]);
                        sumNearWallThickness += thickness[ppPointI]/s;
                    }
                }
            }

            label totNPoints = returnReduce(meshPoints.size(), sumOp<label>());

            // For empty patches, totNPoints is 0.
            scalar avgThickness = 0;
            scalar avgNearWallThickness = 0;

            if (totNPoints > 0)
            {
                //reduce(maxThickness, maxOp<scalar>());
                //reduce(minThickness, minOp<scalar>());
                avgThickness =
                    returnReduce(sumThickness, sumOp<scalar>())
                  / totNPoints;
                avgNearWallThickness =
                    returnReduce(sumNearWallThickness, sumOp<scalar>())
                  / totNPoints;
            }

            Info<< setf(ios_base::left) << setw(maxPatchNameLen)
                << patches[patchI].name() << setprecision(3)
                << " " << setw(8)
                << returnReduce(patches[patchI].size(), sumOp<scalar>())
                << " " << setw(6) << layerParams.numLayers()[patchI]
                << " " << setw(8) << avgNearWallThickness
                << "  " << setw(8) << avgThickness
                //<< " " << setw(8) << minThickness
                //<< " " << setw(8) << maxThickness
                << endl;
        }
        Info<< endl;
    }


    // Calculate wall to medial axis distance for smoothing displacement
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    pointScalarField pointMedialDist
    (
        IOobject
        (
            "pointMedialDist",
            meshRefiner_.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        meshMover().pMesh(),
        dimensionedScalar("pointMedialDist", dimless, 0.0)
    );

    pointVectorField dispVec
    (
        IOobject
        (
            "dispVec",
            meshRefiner_.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        meshMover().pMesh(),
        dimensionedVector("dispVec", dimless, vector::zero)
    );

    pointScalarField medialRatio
    (
        IOobject
        (
            "medialRatio",
            meshRefiner_.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        meshMover().pMesh(),
        dimensionedScalar("medialRatio", dimless, 0.0)
    );

    // Setup information for medial axis smoothing. Calculates medial axis
    // and a smoothed displacement direction.
    // - pointMedialDist : distance to medial axis
    // - dispVec : normalised direction of nearest displacement
    // - medialRatio : ratio of medial distance to wall distance.
    //   (1 at wall, 0 at medial axis)
    medialAxisSmoothingInfo
    (
        meshMover,
        layerParams.nSmoothNormals(),
        layerParams.nSmoothSurfaceNormals(),
        layerParams.minMedianAxisAngleCos(),

        dispVec,
        medialRatio,
        pointMedialDist
    );



    // Saved old points
    pointField oldPoints(mesh.points());

    // Last set of topology changes. (changing mesh clears out polyTopoChange)
    polyTopoChange savedMeshMod(mesh.boundaryMesh().size());

    boolList flaggedCells;
    boolList flaggedFaces;

    for (label iteration = 0; iteration < layerParams.nLayerIter(); iteration++)
    {
        Info<< nl
            << "Layer addition iteration " << iteration << nl
            << "--------------------------" << endl;


        // Unset the extrusion at the pp.
        const dictionary& meshQualityDict =
        (
            iteration < layerParams.nRelaxedIter()
          ? motionDict
          : motionDict.subDict("relaxed")
        );

        if (iteration >= layerParams.nRelaxedIter())
        {
            Info<< "Switched to relaxed meshQuality constraints." << endl;
        }



        // Make sure displacement is equal on both sides of coupled patches.
        syncPatchDisplacement
        (
            meshMover,
            minThickness,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Displacement acc. to pointnormals
        getPatchDisplacement
        (
            meshMover,
            thickness,
            minThickness,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Shrink mesh by displacement value first.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        {
            pointField oldPatchPos(pp().localPoints());

            //// Laplacian displacement shrinking.
            //shrinkMeshDistance
            //(
            //    debug,
            //    meshMover,
            //    -patchDisp,     // Shrink in opposite direction of addedPoints
            //    layerParams.nSmoothDisp(),
            //    layerParams.nSnap()
            //);

            // Medial axis based shrinking
            shrinkMeshMedialDistance
            (
                meshMover(),
                meshQualityDict,

                layerParams.nSmoothThickness(),
                layerParams.maxThicknessToMedialRatio(),
                nAllowableErrors,
                layerParams.nSnap(),
                layerParams.layerTerminationCos(),

                thickness,
                minThickness,

                dispVec,
                medialRatio,
                pointMedialDist,

                extrudeStatus,
                patchDisp,
                patchNLayers
            );

            // Update patchDisp (since not all might have been honoured)
            patchDisp = oldPatchPos - pp().localPoints();
        }

        // Truncate displacements that are too small (this will do internal
        // ones, coupled ones have already been truncated by syncPatch)
        faceSet dummySet(mesh, "wrongPatchFaces", 0);
        truncateDisplacement
        (
            meshMover(),
            minThickness,
            dummySet,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );


        // Dump to .obj file for debugging.
        if (debug)
        {
            dumpDisplacement
            (
                mesh.time().path()/"layer",
                pp(),
                patchDisp,
                extrudeStatus
            );

            const_cast<Time&>(mesh.time())++;
            Info<< "Writing shrunk mesh to " << meshRefiner_.timeName() << endl;

            // See comment in autoSnapDriver why we should not remove meshPhi
            // using mesh.clearPout().

            mesh.write();
        }


        // Mesh topo change engine
        polyTopoChange meshMod(mesh);

        // Grow layer of cells on to patch. Handles zero sized displacement.
        addPatchCellLayer addLayer(mesh);

        // Determine per point/per face number of layers to extrude. Also
        // handles the slow termination of layers when going switching layers

        labelList nPatchPointLayers(pp().nPoints(),-1);
        labelList nPatchFaceLayers(pp().localFaces().size(),-1);
        setupLayerInfoTruncation
        (
            meshMover(),
            patchNLayers,
            extrudeStatus,
            layerParams.nBufferCellsNoExtrude(),
            nPatchPointLayers,
            nPatchFaceLayers
        );

        // Calculate displacement for first layer for addPatchLayer.
        // (first layer = layer of cells next to the original mesh)
        vectorField firstDisp(patchNLayers.size(), vector::zero);

        forAll(patchNLayers, i)
        {
            if (patchNLayers[i] > 0)
            {
                if (expansionRatio[i] == 1.0)
                {
                    firstDisp[i] = patchDisp[i]/nPatchPointLayers[i];
                }
                else
                {
                    label nLay = nPatchPointLayers[i];
                    scalar h =
                        pow(expansionRatio[i], nLay - 1)
                      * (1.0 - expansionRatio[i])
                      / (1.0 - pow(expansionRatio[i], nLay));
                    firstDisp[i] = h*patchDisp[i];
                }
            }
        }

        scalarField invExpansionRatio = 1.0 / expansionRatio;

        // Add topo regardless of whether extrudeStatus is extruderemove.
        // Not add layer if patchDisp is zero.
        addLayer.setRefinement
        (
            invExpansionRatio,
            pp(),
            nPatchFaceLayers,   // layers per face
            nPatchPointLayers,  // layers per point
            firstDisp,          // thickness of layer nearest internal mesh
            meshMod
        );

        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        // Store mesh changes for if mesh is correct.
        savedMeshMod = meshMod;


        // With the stored topo changes we create a new mesh so we can
        // undo if neccesary.

        autoPtr<fvMesh> newMeshPtr;
        autoPtr<mapPolyMesh> map = meshMod.makeMesh
        (
            newMeshPtr,
            IOobject
            (
                //mesh.name()+"_layer",
                mesh.name(),
                static_cast<polyMesh&>(mesh).instance(),
                mesh.time(),  // register with runTime
                static_cast<polyMesh&>(mesh).readOpt(),
                static_cast<polyMesh&>(mesh).writeOpt()
            ),              // io params from original mesh but new name
            mesh,           // original mesh
            true            // parallel sync
        );
        fvMesh& newMesh = newMeshPtr();

        //?neccesary? Update fields
        newMesh.updateMesh(map);

        if (meshRefiner_.overwrite())
        {
            newMesh.setInstance(meshRefiner_.oldInstance());
        }

        // Update numbering on addLayer:
        // - cell/point labels to be newMesh.
        // - patchFaces to remain in oldMesh order.
        addLayer.updateMesh
        (
            map,
            identity(pp().size()),
            identity(pp().nPoints())
        );

        // Collect layer faces and cells for outside loop.
        getLayerCellsFaces
        (
            newMesh,
            addLayer,
            flaggedCells,
            flaggedFaces
        );


        if (debug)
        {
            Info<< "Writing layer mesh to " << meshRefiner_.timeName() << endl;
            newMesh.write();
            cellSet addedCellSet
            (
                newMesh,
                "addedCells",
                findIndices(flaggedCells, true)
            );
            Info<< "Writing "
                << returnReduce(addedCellSet.size(), sumOp<label>())
                << " added cells to cellSet "
                << addedCellSet.name() << endl;
            addedCellSet.write();

            faceSet layerFacesSet
            (
                newMesh,
                "layerFaces",
                findIndices(flaggedCells, true)
            );
            Info<< "Writing "
                << returnReduce(layerFacesSet.size(), sumOp<label>())
                << " faces inside added layer to faceSet "
                << layerFacesSet.name() << endl;
            layerFacesSet.write();
        }


        label nTotChanged = checkAndUnmark
        (
            addLayer,
            meshQualityDict,
            pp(),
            newMesh,

            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        Info<< "Extruding " << countExtrusion(pp, extrudeStatus)
            << " out of " << returnReduce(pp().size(), sumOp<label>())
            << " faces. Removed extrusion at " << nTotChanged << " faces."
            << endl;

        if (nTotChanged == 0)
        {
            break;
        }

        // Reset mesh points and start again
        meshMover().movePoints(oldPoints);
        meshMover().correct();

        Info<< endl;
    }


    // At this point we have a (shrunk) mesh and a set of topology changes
    // which will make a valid mesh with layer. Apply these changes to the
    // current mesh.

    // Apply the stored topo changes to the current mesh.
    autoPtr<mapPolyMesh> map = savedMeshMod.changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh.clearOut();
    }

    if (meshRefiner_.overwrite())
    {
        mesh.setInstance(meshRefiner_.oldInstance());
    }

    meshRefiner_.updateMesh(map, labelList(0));


    // Do final balancing
    // ~~~~~~~~~~~~~~~~~~

    if (Pstream::parRun())
    {
        Info<< nl
            << "Doing final balancing" << nl
            << "---------------------" << nl
            << endl;

        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        // Balance. No restriction on face zones and baffles.
        autoPtr<mapDistributePolyMesh> map = meshRefiner_.balance
        (
            false,
            false,
            scalarField(mesh.nCells(), 1.0),
            decomposer,
            distributor
        );

        // Re-distribute flag of layer faces and cells
        map().distributeCellData(flaggedCells);
        map().distributeFaceData(flaggedFaces);
    }


    // Write mesh
    // ~~~~~~~~~~

    cellSet addedCellSet(mesh, "addedCells", findIndices(flaggedCells, true));
    Info<< "Writing "
        << returnReduce(addedCellSet.size(), sumOp<label>())
        << " added cells to cellSet "
        << addedCellSet.name() << endl;
    addedCellSet.write();

    faceSet layerFacesSet(mesh, "layerFaces", findIndices(flaggedFaces, true));
    Info<< "Writing "
        << returnReduce(layerFacesSet.size(), sumOp<label>())
        << " faces inside added layer to faceSet "
        << layerFacesSet.name() << endl;
    layerFacesSet.write();
}


void Foam::autoLayerDriver::doLayers
(
    const dictionary& shrinkDict,
    const dictionary& motionDict,
    const layerParameters& layerParams,
    const bool preBalance,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl
        << "Shrinking and layer addition phase" << nl
        << "----------------------------------" << nl
        << endl;

    Info<< "Using mesh parameters " << motionDict << nl << endl;

    // Merge coplanar boundary faces
    mergePatchFacesUndo(layerParams, motionDict);

    // Per patch the number of layers (0 if no layer)
    const labelList& numLayers = layerParams.numLayers();

    // Patches that need to get a layer
    DynamicList<label> patchIDs(numLayers.size());
    label nFacesWithLayers = 0;
    forAll(numLayers, patchI)
    {
        if (numLayers[patchI] > 0)
        {
            patchIDs.append(patchI);
            nFacesWithLayers += mesh.boundaryMesh()[patchI].size();
        }
    }
    patchIDs.shrink();

    if (returnReduce(nFacesWithLayers, sumOp<label>()) == 0)
    {
        Info<< nl << "No layers to generate ..." << endl;
    }
    else
    {
        // Check that outside of mesh is not multiply connected.
        checkMeshManifold();

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


        // Balance
        if (Pstream::parRun() && preBalance)
        {
            Info<< nl
                << "Doing initial balancing" << nl
                << "-----------------------" << nl
                << endl;

            scalarField cellWeights(mesh.nCells(), 1);
            forAll(numLayers, patchI)
            {
                if (numLayers[patchI] > 0)
                {
                    const polyPatch& pp = mesh.boundaryMesh()[patchI];
                    forAll(pp.faceCells(), i)
                    {
                        cellWeights[pp.faceCells()[i]] += numLayers[patchI];
                    }
                }
            }

            // Balance mesh (and meshRefinement). No restriction on face zones
            // and baffles.
            autoPtr<mapDistributePolyMesh> map = meshRefiner_.balance
            (
                false,
                false,
                cellWeights,
                decomposer,
                distributor
            );

            //{
            //    globalIndex globalCells(mesh.nCells());
            //
            //    Info<< "** Distribution after balancing:" << endl;
            //    for (label procI = 0; procI < Pstream::nProcs(); procI++)
            //    {
            //        Info<< "    " << procI << '\t'
            //            << globalCells.localSize(procI) << endl;
            //    }
            //    Info<< endl;
            //}
        }


        // Do all topo changes
        addLayers
        (
            layerParams,
            motionDict,
            patchIDs,
            nInitErrors,
            decomposer,
            distributor
        );
    }
}


// ************************************************************************* //
