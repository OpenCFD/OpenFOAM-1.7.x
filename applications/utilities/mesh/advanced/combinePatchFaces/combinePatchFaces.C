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
    Checks for multiple patch faces on same cell and combines them. These
    result from e.g. refined neighbouring cells getting removed, leaving 4
    exposed faces with same owner.

    Rules for merging:
    - only boundary faces (since multiple internal faces between two cells
      not allowed anyway)
    - faces have to have same owner
    - faces have to be connected via edge which are not features (so angle
      between them < feature angle)
    - outside of faces has to be single loop
    - outside of face should not be (or just slightly) concave (so angle
      between consecutive edges < concaveangle

    E.g. to allow all faces on same patch to be merged:

        combinePatchFaces .. cavity 180 -concaveAngle 90

\*---------------------------------------------------------------------------*/

#include "PstreamReduceOps.H"
#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "polyModifyFace.H"
#include "polyAddFace.H"
#include "combineFaces.H"
#include "removePoints.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "mathematicalConstants.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Sin of angle between two consecutive edges on a face. If sin(angle) larger
// than this the face will be considered concave.
const scalar defaultConcaveAngle = 30;


// Same check as snapMesh
void checkSnapMesh
(
    const Time& runTime,
    const polyMesh& mesh,
    labelHashSet& wrongFaces
)
{
    IOdictionary snapDict
    (
        IOobject
        (
            "snapMeshDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Max nonorthogonality allowed
    scalar maxNonOrtho(readScalar(snapDict.lookup("maxNonOrtho")));
    // Max concaveness allowed.
    scalar maxConcave(readScalar(snapDict.lookup("maxConcave")));
    // Min volume allowed (factor of minimum cellVolume)
    scalar relMinVol(readScalar(snapDict.lookup("minVol")));
    const scalar minCellVol = min(mesh.cellVolumes());
    const scalar minPyrVol = relMinVol*minCellVol;
    // Min area
    scalar minArea(readScalar(snapDict.lookup("minArea")));

    if (maxNonOrtho < 180.0-SMALL)
    {
        Pout<< "Checking non orthogonality" << endl;

        label nOldSize = wrongFaces.size();
        mesh.setNonOrthThreshold(maxNonOrtho);
        mesh.checkFaceOrthogonality(false, &wrongFaces);

        Pout<< "Detected " << wrongFaces.size() - nOldSize
            << " faces with non-orthogonality > " << maxNonOrtho << " degrees"
            << endl;
    }

    if (minPyrVol > -GREAT)
    {
        Pout<< "Checking face pyramids" << endl;

        label nOldSize = wrongFaces.size();
        mesh.checkFacePyramids(false, minPyrVol, &wrongFaces);
        Pout<< "Detected additional " << wrongFaces.size() - nOldSize
            << " faces with illegal face pyramids" << endl;
    }

    if (maxConcave < 180.0-SMALL)
    {
        Pout<< "Checking face angles" << endl;

        label nOldSize = wrongFaces.size();
        mesh.checkFaceAngles(false, maxConcave, &wrongFaces);
        Pout<< "Detected additional " << wrongFaces.size() - nOldSize
            << " faces with concavity > " << maxConcave << " degrees"
            << endl;
    }

    if (minArea > -SMALL)
    {
        Pout<< "Checking face areas" << endl;

        label nOldSize = wrongFaces.size();

        const scalarField magFaceAreas = mag(mesh.faceAreas());

        forAll(magFaceAreas, faceI)
        {
            if (magFaceAreas[faceI] < minArea)
            {
                wrongFaces.insert(faceI);
            }
        }
        Pout<< "Detected additional " << wrongFaces.size() - nOldSize
            << " faces with area < " << minArea << " m^2" << endl;
    }
}


// Merge faces on the same patch (usually from exposing refinement)
// Can undo merges if these cause problems.
label mergePatchFaces
(
    const scalar minCos,
    const scalar concaveSin,
    const bool snapMeshDict,
    const Time& runTime,
    polyMesh& mesh
)
{
    // Patch face merging engine
    combineFaces faceCombiner(mesh);

    // Get all sets of faces that can be merged
    labelListList allFaceSets(faceCombiner.getMergeSets(minCos, concaveSin));

    label nFaceSets = returnReduce(allFaceSets.size(), sumOp<label>());

    Info<< "Merging " << nFaceSets << " sets of faces." << endl;

    if (nFaceSets > 0)
    {
        // Store the faces of the face sets
        List<faceList> allFaceSetsFaces(allFaceSets.size());
        forAll(allFaceSets, setI)
        {
            allFaceSetsFaces[setI] = UIndirectList<face>
            (
                mesh.faces(),
                allFaceSets[setI]
            );
        }

        autoPtr<mapPolyMesh> map;
        {
            // Topology changes container
            polyTopoChange meshMod(mesh);

            // Merge all faces of a set into the first face of the set.
            faceCombiner.setRefinement(allFaceSets, meshMod);

            // Change the mesh (no inflation)
            map = meshMod.changeMesh(mesh, false, true);

            // Update fields
            mesh.updateMesh(map);

            // Move mesh (since morphing does not do this)
            if (map().hasMotionPoints())
            {
                mesh.movePoints(map().preMotionPoints());
            }
            else
            {
                // Delete mesh volumes. No other way to do this?
                mesh.clearOut();
            }
        }


        // Check for errors and undo
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        // Faces in error.
        labelHashSet errorFaces;

        if (snapMeshDict)
        {
            checkSnapMesh(runTime, mesh, errorFaces);
        }
        else
        {
            mesh.checkFacePyramids(false, -SMALL, &errorFaces);
        }

        // Sets where the master is in error
        labelHashSet errorSets;

        forAll(allFaceSets, setI)
        {
            label newMasterI = map().reverseFaceMap()[allFaceSets[setI][0]];

            if (errorFaces.found(newMasterI))
            {
                errorSets.insert(setI);
            }
        }
        label nErrorSets = returnReduce(errorSets.size(), sumOp<label>());

        Info<< "Detected " << nErrorSets
            << " error faces on boundaries that have been merged."
            << " These will be restored to their original faces."
            << endl;

        if (nErrorSets > 0)
        {
            // Renumber stored faces to new vertex numbering.
            forAllConstIter(labelHashSet, errorSets, iter)
            {
                label setI = iter.key();

                faceList& setFaceVerts = allFaceSetsFaces[setI];

                forAll(setFaceVerts, i)
                {
                    inplaceRenumber(map().reversePointMap(), setFaceVerts[i]);

                    // Debug: check that all points are still there.
                    forAll(setFaceVerts[i], j)
                    {
                        label newVertI = setFaceVerts[i][j];

                        if (newVertI < 0)
                        {
                            FatalErrorIn("mergePatchFaces")
                                << "In set:" << setI << " old face labels:"
                                << allFaceSets[setI] << " new face vertices:"
                                << setFaceVerts[i] << " are unmapped vertices!"
                                << abort(FatalError);
                        }
                    }
                }
            }


            // Topology changes container
            polyTopoChange meshMod(mesh);


            // Restore faces
            forAllConstIter(labelHashSet, errorSets, iter)
            {
                label setI = iter.key();

                const labelList& setFaces = allFaceSets[setI];
                const faceList& setFaceVerts = allFaceSetsFaces[setI];

                label newMasterI = map().reverseFaceMap()[setFaces[0]];

                // Restore. Get face properties.

                label own = mesh.faceOwner()[newMasterI];
                label zoneID = mesh.faceZones().whichZone(newMasterI);
                bool zoneFlip = false;
                if (zoneID >= 0)
                {
                    const faceZone& fZone = mesh.faceZones()[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(newMasterI)];
                }
                label patchI = mesh.boundaryMesh().whichPatch(newMasterI);


                Pout<< "Restoring new master face " << newMasterI
                    << " to vertices " << setFaceVerts[0] << endl;

                // Modify the master face.
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        setFaceVerts[0],                // original face
                        newMasterI,                     // label of face
                        own,                            // owner
                        -1,                             // neighbour
                        false,                          // face flip
                        patchI,                         // patch for face
                        false,                          // remove from zone
                        zoneID,                         // zone for face
                        zoneFlip                        // face flip in zone
                    )
                );


                // Add the previously removed faces
                for (label i = 1; i < setFaces.size(); i++)
                {
                    Pout<< "Restoring removed face " << setFaces[i]
                        << " with vertices " << setFaceVerts[i] << endl;

                    meshMod.setAction
                    (
                        polyAddFace
                        (
                            setFaceVerts[i],        // vertices
                            own,                    // owner,
                            -1,                     // neighbour,
                            -1,                     // masterPointID,
                            -1,                     // masterEdgeID,
                            newMasterI,             // masterFaceID,
                            false,                  // flipFaceFlux,
                            patchI,                 // patchID,
                            zoneID,                 // zoneID,
                            zoneFlip                // zoneFlip
                        )
                    );
                }
            }

            // Change the mesh (no inflation)
            map = meshMod.changeMesh(mesh, false, true);

            // Update fields
            mesh.updateMesh(map);

            // Move mesh (since morphing does not do this)
            if (map().hasMotionPoints())
            {
                mesh.movePoints(map().preMotionPoints());
            }
            else
            {
                // Delete mesh volumes. No other way to do this?
                mesh.clearOut();
            }
        }
    }
    else
    {
        Info<< "No faces merged ..." << endl;
    }

    return nFaceSets;
}


// Remove points not used by any face or points used by only two faces where
// the edges are in line
label mergeEdges(const scalar minCos, polyMesh& mesh)
{
    Info<< "Merging all points on surface that" << nl
        << "- are used by only two boundary faces and" << nl
        << "- make an angle with a cosine of more than " << minCos
        << "." << nl << endl;

    // Point removal analysis engine
    removePoints pointRemover(mesh);

    // Count usage of points
    boolList pointCanBeDeleted;
    label nRemove = pointRemover.countPointUsage(minCos, pointCanBeDeleted);

    if (nRemove > 0)
    {
        Info<< "Removing " << nRemove
            << " straight edge points ..." << endl;

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
            // Delete mesh volumes. No other way to do this?
            mesh.clearOut();
        }
    }
    else
    {
        Info<< "No straight edges simplified and no points removed ..." << endl;
    }

    return nRemove;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("feature angle [0..180]");
    argList::validOptions.insert("concaveAngle", "[0..180]");
    argList::validOptions.insert("snapMesh", "");
    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
    runTime.functionObjects().off();
#   include "createPolyMesh.H"
    const word oldInstance = mesh.pointsInstance();

    scalar featureAngle(readScalar(IStringStream(args.additionalArgs()[0])()));

    scalar minCos = Foam::cos(featureAngle*mathematicalConstant::pi/180.0);

    scalar concaveAngle = defaultConcaveAngle;
    args.optionReadIfPresent("concaveAngle", concaveAngle);

    scalar concaveSin = Foam::sin(concaveAngle*mathematicalConstant::pi/180.0);

    bool snapMeshDict = args.optionFound("snapMesh");
    bool overwrite = args.optionFound("overwrite");

    Info<< "Merging all faces of a cell" << nl
        << "    - which are on the same patch" << nl
        << "    - which make an angle < " << featureAngle << " degrees"
        << nl
        << "      (cos:" << minCos << ')' << nl
        << "    - even when resulting face becomes concave by more than "
        << concaveAngle << " degrees" << nl
        << "      (sin:" << concaveSin << ')' << nl
        << endl;

    if (!overwrite)
    {
        runTime++;
    }

    // Merge faces on same patch
    label nChanged = mergePatchFaces
    (
        minCos,
        concaveSin,
        snapMeshDict,
        runTime,
        mesh
    );

    // Merge points on straight edges and remove unused points
    if (snapMeshDict)
    {
        Info<< "Merging all 'loose' points on surface edges"
            << ", regardless of the angle they make." << endl;

        // Surface bnound to be used to extrude. Merge all loose points.
        nChanged += mergeEdges(-1, mesh);
    }
    else
    {
        nChanged += mergeEdges(minCos, mesh);
    }

    if (nChanged > 0)
    {
        if (overwrite)
        {
            mesh.setInstance(oldInstance);
        }

        Info<< "Writing morphed mesh to time " << runTime.timeName() << endl;

        mesh.write();
    }
    else
    {
        Info<< "Mesh unchanged." << endl;
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
