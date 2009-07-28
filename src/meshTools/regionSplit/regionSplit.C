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

\*---------------------------------------------------------------------------*/

#include "regionSplit.H"
#include "cyclicPolyPatch.H"
#include "processorPolyPatch.H"
#include "globalIndex.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(regionSplit, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Handle (non-processor) coupled faces.
void Foam::regionSplit::transferCoupledFaceRegion
(
    const label faceI,
    const label otherFaceI,

    labelList& faceRegion,
    DynamicList<label>& newChangedFaces
) const
{
    if (faceRegion[faceI] >= 0)
    {
        if (faceRegion[otherFaceI] == -1)
        {
            faceRegion[otherFaceI] = faceRegion[faceI];
            newChangedFaces.append(otherFaceI);
        }
        else if (faceRegion[otherFaceI] == -2)
        {
            // otherFaceI blocked but faceI is not. Is illegal for coupled
            // faces, not for explicit connections.
        }
        else if (faceRegion[otherFaceI] != faceRegion[faceI])
        {
            FatalErrorIn
            (
                  "regionSplit::transferCoupledFaceRegion"
                  "(const label, const label, labelList&, labelList&) const"
              )   << "Problem : coupled face " << faceI
                  << " on patch " << mesh_.boundaryMesh().whichPatch(faceI)
                  << " has region " << faceRegion[faceI]
                  << " but coupled face " << otherFaceI
                  << " has region " << faceRegion[otherFaceI]
                  << endl
                  << "Is your blocked faces specification"
                  << " synchronized across coupled boundaries?"
                  << abort(FatalError);
        }
    }
    else if (faceRegion[faceI] == -1)
    {
        if (faceRegion[otherFaceI] >= 0)
        {
            faceRegion[faceI] = faceRegion[otherFaceI];
            newChangedFaces.append(faceI);
        }
        else if (faceRegion[otherFaceI] == -2)
        {
            // otherFaceI blocked but faceI is not. Is illegal for coupled
            // faces, not for explicit connections.
        }
    }
}


void Foam::regionSplit::fillSeedMask
(
    const List<labelPair>& explicitConnections,
    labelList& cellRegion,
    labelList& faceRegion,
    const label seedCellID,
    const label markValue
) const
{
    // Do seed cell
    cellRegion[seedCellID] = markValue;


    // Collect faces on seed cell
    const cell& cFaces = mesh_.cells()[seedCellID];

    label nFaces = 0;

    labelList changedFaces(cFaces.size());

    forAll(cFaces, i)
    {
        label faceI = cFaces[i];

        if (faceRegion[faceI] == -1)
        {
            faceRegion[faceI] = markValue;
            changedFaces[nFaces++] = faceI;
        }
    }
    changedFaces.setSize(nFaces);


    // Loop over changed faces. MeshWave in small.

    while (changedFaces.size())
    {
        //if (debug)
        //{
        //    Pout<< "regionSplit::fillSeedMask : changedFaces:"
        //        << changedFaces.size() << endl;
        //}

        DynamicList<label> changedCells(changedFaces.size());

        forAll(changedFaces, i)
        {
            label faceI = changedFaces[i];

            label own = mesh_.faceOwner()[faceI];

            if (cellRegion[own] == -1)
            {
                cellRegion[own] = markValue;
                changedCells.append(own);
            }

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];

                if (cellRegion[nei] == -1)
                {
                    cellRegion[nei] = markValue;
                    changedCells.append(nei);
                }
            }
        }


        //if (debug)
        //{
        //    Pout<< "regionSplit::fillSeedMask : changedCells:"
        //        << changedCells.size() << endl;
        //}

        // Loop over changedCells and collect faces
        DynamicList<label> newChangedFaces(changedCells.size());

        forAll(changedCells, i)
        {
            label cellI = changedCells[i];

            const cell& cFaces = mesh_.cells()[cellI];

            forAll(cFaces, cFaceI)
            {
                label faceI = cFaces[cFaceI];

                if (faceRegion[faceI] == -1)
                {
                    faceRegion[faceI] = markValue;
                    newChangedFaces.append(faceI);
                }
            }
        }


        //if (debug)
        //{
        //    Pout<< "regionSplit::fillSeedMask : changedFaces before sync:"
        //        << changedFaces.size() << endl;
        //}


        // Check for changes to any locally coupled face.
        // Global connections are done later.

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (isA<cyclicPolyPatch>(pp))
            {
                label faceI = pp.start();

                label halfSz = pp.size()/2;

                for (label i = 0; i < halfSz; i++)
                {
                    label otherFaceI = refCast<const cyclicPolyPatch>(pp)
                        .transformGlobalFace(faceI);

                    transferCoupledFaceRegion
                    (
                        faceI,
                        otherFaceI,
                        faceRegion,
                        newChangedFaces
                    );

                    faceI++;
                }
            }
        }
        forAll(explicitConnections, i)
        {
            transferCoupledFaceRegion
            (
                explicitConnections[i][0],
                explicitConnections[i][1],
                faceRegion,
                newChangedFaces
            );
        }

        //if (debug)
        //{
        //    Pout<< "regionSplit::fillSeedMask : changedFaces after sync:"
        //        << newChangedFaces.size() << endl;
        //}

        changedFaces.transfer(newChangedFaces);
    }
}


Foam::label Foam::regionSplit::calcRegionSplit
(
    const boolList& blockedFace,
    const List<labelPair>& explicitConnections,

    labelList& cellRegion
) const
{
    if (debug)
    {
        if (blockedFace.size())
        {
            // Check that blockedFace is synced.
            boolList syncBlockedFace(blockedFace);
            syncTools::swapFaceList(mesh_, syncBlockedFace, false);

            forAll(syncBlockedFace, faceI)
            {
                if (syncBlockedFace[faceI] != blockedFace[faceI])
                {
                    FatalErrorIn
                    (
                        "regionSplit::calcRegionSplit(..)"
                    )   << "Face " << faceI << " not synchronised. My value:"
                        << blockedFace[faceI] << "  coupled value:"
                        << syncBlockedFace[faceI]
                        << abort(FatalError);
                }
            }
        }
    }

    // Region per face.
    // -1 unassigned
    // -2 blocked
    labelList faceRegion(mesh_.nFaces(), -1);

    if (blockedFace.size())
    {
        forAll(blockedFace, faceI)
        {
            if (blockedFace[faceI])
            {
                faceRegion[faceI] = -2;
            }
        }
    }


    // Assign local regions
    // ~~~~~~~~~~~~~~~~~~~~

    // Start with region 0
    label nRegions = 0;

    label unsetCellI = 0;

    do
    {
        // Find first unset cell

        for (; unsetCellI < mesh_.nCells(); unsetCellI++)
        {
            if (cellRegion[unsetCellI] == -1)
            {
                break;
            }
        }

        if (unsetCellI >= mesh_.nCells())
        {
            break;
        }

        fillSeedMask
        (
            explicitConnections,
            cellRegion,
            faceRegion,
            unsetCellI,
            nRegions
        );

        // Current unsetCell has now been handled. Go to next region.
        nRegions++;
        unsetCellI++;
    }
    while(true);


    if (debug)
    {
        forAll(cellRegion, cellI)
        {
            if (cellRegion[cellI] < 0)
            {
                FatalErrorIn("regionSplit::calcRegionSplit(..)")
                    << "cell:" << cellI << " region:" << cellRegion[cellI]
                    << abort(FatalError);
            }
        }

        forAll(faceRegion, faceI)
        {
            if (faceRegion[faceI] == -1)
            {
                FatalErrorIn("regionSplit::calcRegionSplit(..)")
                    << "face:" << faceI << " region:" << faceRegion[faceI]
                    << abort(FatalError);
            }
        }
    }



    // Assign global regions
    // ~~~~~~~~~~~~~~~~~~~~~
    // Offset local regions to create unique global regions.

    globalIndex globalRegions(nRegions);


    // Merge global regions
    // ~~~~~~~~~~~~~~~~~~~~
    // Regions across non-blocked proc patches get merged.
    // This will set merged global regions to be the min of both.
    // (this will create gaps in the global region list so they will get
    // merged later on)

    // Map from global to merged global
    labelList mergedGlobal(identity(globalRegions.size()));


    // See if any regions get merged. Only nessecary for parallel
    while (Pstream::parRun())
    {
        if (debug)
        {
            Pout<< nl << "-- Starting Iteration --" << endl;
        }

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        // Send global regions across (or -2 if blocked face)
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (isA<processorPolyPatch>(pp))
            {
                labelList myGlobalRegions(pp.size());

                label faceI = pp.start();

                forAll(pp, i)
                {
                    if (faceRegion[faceI] < 0)
                    {
                        myGlobalRegions[i] = faceRegion[faceI];
                    }
                    else
                    {
                        myGlobalRegions[i] = mergedGlobal
                        [globalRegions.toGlobal(faceRegion[faceI])];
                    }

                    faceI++;
                }

                OPstream toProcNbr
                (
                    Pstream::blocking,
                    refCast<const processorPolyPatch>(pp).neighbProcNo()
                );

                toProcNbr << myGlobalRegions;
            }
        }


        // Receive global regions

        label nMerged = 0;

        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (isA<processorPolyPatch>(pp))
            {
                const processorPolyPatch& procPp =
                    refCast<const processorPolyPatch>(pp);

                IPstream fromProcNbr(Pstream::blocking, procPp.neighbProcNo());

                labelList nbrRegions(fromProcNbr);


                // Compare with my regions to see which get merged.

                label faceI = pp.start();

                forAll(pp, i)
                {
                    if
                    (
                        faceRegion[faceI] < 0
                     || nbrRegions[i] < 0
                    )
                    {
                        if (faceRegion[faceI] != nbrRegions[i])
                        {
                            FatalErrorIn("regionSplit::calcRegionSplit(..)")
                                << "On patch:" << pp.name()
                                << " face:" << faceI
                                << " my local region:" << faceRegion[faceI]
                                << " neighbouring region:"
                                << nbrRegions[i] << nl
                                << "Maybe your blockedFaces are not"
                                << " synchronized across coupled faces?"
                                << abort(FatalError);
                        }
                    }
                    else
                    {
                        label uncompactGlobal =
                            globalRegions.toGlobal(faceRegion[faceI]);

                        label myGlobal = mergedGlobal[uncompactGlobal];

                        if (myGlobal != nbrRegions[i])
                        {
                            label minRegion = min(myGlobal, nbrRegions[i]);

                            if (debug)
                            {
                                Pout<< "Merging region " << myGlobal
                                    << " (on proc " << Pstream::myProcNo()
                                    << ") and region " << nbrRegions[i]
                                    << " (on proc " << procPp.neighbProcNo()
                                    << ") into region " << minRegion << endl;
                            }

                            mergedGlobal[uncompactGlobal] = minRegion;
                            mergedGlobal[myGlobal] = minRegion;
                            mergedGlobal[nbrRegions[i]] = minRegion;

                            nMerged++;
                        }
                    }

                    faceI++;
                }
            }
        }


        reduce(nMerged, sumOp<label>());

        if (debug)
        {
            Pout<< "nMerged:" << nMerged << endl;
        }

        if (nMerged == 0)
        {
            break;
        }

        // Merge the compacted regions.
        Pstream::listCombineGather(mergedGlobal, minEqOp<label>());
        Pstream::listCombineScatter(mergedGlobal);
    }


    // Compact global regions
    // ~~~~~~~~~~~~~~~~~~~~~~

    // All procs will have the same global mergedGlobal region.
    // There might be gaps in it however so compact.

    labelList mergedToCompacted(globalRegions.size(), -1);

    label compactI = 0;

    forAll(mergedGlobal, i)
    {
        label merged = mergedGlobal[i];

        if (mergedToCompacted[merged] == -1)
        {
            mergedToCompacted[merged] = compactI++;
        }
    }

    if (debug)
    {
        Pout<< "Compacted down to " << compactI << " regions." << endl;
    }

    // Renumber cellRegion to be global regions
    forAll(cellRegion, cellI)
    {
        label region = cellRegion[cellI];

        if (region >= 0)
        {
            label merged = mergedGlobal[globalRegions.toGlobal(region)];

            cellRegion[cellI] = mergedToCompacted[merged];
        }
    }

    return compactI;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSplit::regionSplit(const polyMesh& mesh)
:
    labelList(mesh.nCells(), -1),
    mesh_(mesh),
    nRegions_(calcRegionSplit(boolList(0, false), List<labelPair>(0), *this))
{}


Foam::regionSplit::regionSplit
(
    const polyMesh& mesh,
    const boolList& blockedFace
)
:
    labelList(mesh.nCells(), -1),
    mesh_(mesh),
    nRegions_(calcRegionSplit(blockedFace, List<labelPair>(0), *this))
{}


Foam::regionSplit::regionSplit
(
    const polyMesh& mesh,
    const boolList& blockedFace,
    const List<labelPair>& explicitConnections
)
:
    labelList(mesh.nCells(), -1),
    mesh_(mesh),
    nRegions_(calcRegionSplit(blockedFace, explicitConnections, *this))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
