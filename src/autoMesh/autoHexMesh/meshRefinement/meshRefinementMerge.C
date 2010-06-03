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

\*----------------------------------------------------------------------------*/

#include "meshRefinement.H"
#include "combineFaces.H"
#include "polyTopoChange.H"
#include "removePoints.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Merge faces that are in-line.
Foam::label Foam::meshRefinement::mergePatchFaces
(
    const scalar minCos,
    const scalar concaveCos,
    const labelList& patchIDs
)
{
    // Patch face merging engine
    combineFaces faceCombiner(mesh_);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Pick up all candidate cells on boundary
    labelHashSet boundaryCells(mesh_.nFaces()-mesh_.nInternalFaces());

    forAll(patchIDs, i)
    {
        label patchI = patchIDs[i];

        const polyPatch& patch = patches[patchI];

        if (!patch.coupled())
        {
            forAll(patch, i)
            {
                boundaryCells.insert(mesh_.faceOwner()[patch.start()+i]);
            }
        }
    }

    // Get all sets of faces that can be merged
    labelListList mergeSets
    (
        faceCombiner.getMergeSets
        (
            minCos,
            concaveCos,
            boundaryCells
        )
    );

    label nFaceSets = returnReduce(mergeSets.size(), sumOp<label>());

    Info<< "mergePatchFaces : Merging " << nFaceSets
        << " sets of faces." << endl;

    if (nFaceSets > 0)
    {
        // Topology changes container
        polyTopoChange meshMod(mesh_);

        // Merge all faces of a set into the first face of the set. Remove
        // unused points.
        faceCombiner.setRefinement(mergeSets, meshMod);

        // Change the mesh (no inflation)
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh (since morphing does not do this)
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes. No other way to do this?
            mesh_.clearOut();
        }

        if (overwrite())
        {
            mesh_.setInstance(oldInstance());
        }

        faceCombiner.updateMesh(map);

        // Get the kept faces that need to be recalculated.
        // Merging two boundary faces might shift the cell centre
        // (unless the faces are absolutely planar)
        labelHashSet retestFaces(6*mergeSets.size());

        forAll(mergeSets, setI)
        {
            label oldMasterI = mergeSets[setI][0];

            label faceI = map().reverseFaceMap()[oldMasterI];

            // faceI is always uncoupled boundary face
            const cell& cFaces = mesh_.cells()[mesh_.faceOwner()[faceI]];

            forAll(cFaces, i)
            {
                retestFaces.insert(cFaces[i]);
            }
        }
        updateMesh(map, retestFaces.toc());
    }


    return nFaceSets;
}


// Remove points not used by any face or points used by only two faces where
// the edges are in line
Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::mergeEdges
(
    const scalar minCos
)
{
    // Point removal analysis engine
    removePoints pointRemover(mesh_);

    // Count usage of points
    boolList pointCanBeDeleted;
    label nRemove = pointRemover.countPointUsage(minCos, pointCanBeDeleted);

    Info<< "Removing " << nRemove
        << " straight edge points." << endl;

    autoPtr<mapPolyMesh> map;

    if (nRemove > 0)
    {
        // Save my local faces that will change. These changed faces might
        // cause a shift in the cell centre which needs to be retested.
        // Have to do this before changing mesh since point will be removed.
        labelHashSet retestOldFaces(nRemove / Pstream::nProcs());

        {
            const faceList& faces = mesh_.faces();

            forAll(faces, faceI)
            {
                const face& f = faces[faceI];

                forAll(f, fp)
                {
                    if (pointCanBeDeleted[f[fp]])
                    {
                        retestOldFaces.insert(faceI);
                        break;
                    }
                }
            }
        }

        // Topology changes container
        polyTopoChange meshMod(mesh_);

        pointRemover.setRefinement(pointCanBeDeleted, meshMod);

        // Change the mesh (no inflation)
        map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh (since morphing does not do this)
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes. No other way to do this?
            mesh_.clearOut();
        }

        if (overwrite())
        {
            mesh_.setInstance(oldInstance());
        }

        pointRemover.updateMesh(map);

        // Get the kept faces that need to be recalculated.
        labelHashSet retestFaces(6*retestOldFaces.size());

        const cellList& cells = mesh_.cells();

        forAllConstIter(labelHashSet, retestOldFaces, iter)
        {
            label faceI = map().reverseFaceMap()[iter.key()];

            const cell& ownFaces = cells[mesh_.faceOwner()[faceI]];

            forAll(ownFaces, i)
            {
                retestFaces.insert(ownFaces[i]);
            }

            if (mesh_.isInternalFace(faceI))
            {
                const cell& neiFaces = cells[mesh_.faceNeighbour()[faceI]];

                forAll(neiFaces, i)
                {
                    retestFaces.insert(neiFaces[i]);
                }
            }
        }
        updateMesh(map, retestFaces.toc());
    }

    return map;
}


// ************************************************************************* //
