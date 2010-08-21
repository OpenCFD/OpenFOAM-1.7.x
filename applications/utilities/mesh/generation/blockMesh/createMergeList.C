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

#include "blockMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::blockMesh::createMergeList()
{
    Info<< nl << "Creating merge list " << flush;

    labelList MergeList(nPoints_, -1);

    blockMesh& blocks = *this;

    const pointField& blockPoints = topology().points();
    const cellList& blockCells = topology().cells();
    const faceList& blockFaces = topology().faces();
    const labelList& faceOwnerBlocks = topology().faceOwner();

    // For efficiency, create merge pairs in the first pass
    labelListListList glueMergePairs(blockFaces.size());

    const labelList& faceNeighbourBlocks = topology().faceNeighbour();

    forAll(blockFaces, blockFaceLabel)
    {
        label blockPlabel = faceOwnerBlocks[blockFaceLabel];
        const pointField& blockPpoints = blocks[blockPlabel].points();
        const labelList& blockPfaces = blockCells[blockPlabel];

        bool foundFace = false;
        label blockPfaceLabel;
        for
        (
            blockPfaceLabel = 0;
            blockPfaceLabel < blockPfaces.size();
            blockPfaceLabel++
        )
        {
            if
            (
                blockFaces[blockPfaces[blockPfaceLabel]]
             == blockFaces[blockFaceLabel]
            )
            {
                foundFace = true;
                break;
            }
        }

        if (!foundFace)
        {
            FatalErrorIn("blockMesh::createMergeList()")
                << "Cannot find merge face for block " << blockPlabel
                << exit(FatalError);
        };

        const labelListList& blockPfaceFaces =
            blocks[blockPlabel].boundaryPatches()[blockPfaceLabel];

        labelListList& curPairs = glueMergePairs[blockFaceLabel];
        curPairs.setSize(blockPfaceFaces.size());

        // Calculate sqr of the merge tolerance as 1/10th of the min sqr
        // point to point distance on the block face.
        // At the same time merge collated points on the block's faces
        // (removes boundary poles etc.)
        // Co-located points detected by initally taking a constant factor of
        // the size of the block face:
        boundBox bb(blockFaces[blockFaceLabel].points(blockPoints));
        const scalar mergeSqrDist =
            SMALL*magSqr(bb.span())/blockPfaceFaces.size();
        scalar sqrMergeTol = GREAT;

        // This is an N^2 algorithm
        forAll(blockPfaceFaces, blockPfaceFaceLabel)
        {
            const labelList& blockPfaceFacePoints
                = blockPfaceFaces[blockPfaceFaceLabel];

            forAll(blockPfaceFacePoints, blockPfaceFacePointLabel)
            {
                forAll(blockPfaceFacePoints, blockPfaceFacePointLabel2)
                {
                    if (blockPfaceFacePointLabel != blockPfaceFacePointLabel2)
                    {
                        scalar magSqrDist = magSqr
                        (
                            blockPpoints[blockPfaceFacePoints
                            [blockPfaceFacePointLabel]]
                            - blockPpoints[blockPfaceFacePoints
                            [blockPfaceFacePointLabel2]]
                        );

                        if (magSqrDist < mergeSqrDist)
                        {
                            label PpointLabel =
                                blockPfaceFacePoints[blockPfaceFacePointLabel]
                              + blockOffsets_[blockPlabel];

                            label PpointLabel2 =
                                blockPfaceFacePoints[blockPfaceFacePointLabel2]
                              + blockOffsets_[blockPlabel];

                            label minPP2 = min(PpointLabel, PpointLabel2);

                            if (MergeList[PpointLabel] != -1)
                            {
                                minPP2 = min(minPP2, MergeList[PpointLabel]);
                            }

                            if (MergeList[PpointLabel2] != -1)
                            {
                                minPP2 = min(minPP2, MergeList[PpointLabel2]);
                            }

                            MergeList[PpointLabel] = MergeList[PpointLabel2]
                                = minPP2;
                        }
                        else
                        {
                            sqrMergeTol = min(sqrMergeTol, magSqrDist);
                        }
                    }
                }
            }
        }

        sqrMergeTol /= 10.0;


        if (topology().isInternalFace(blockFaceLabel))
        {
            label blockNlabel = faceNeighbourBlocks[blockFaceLabel];
            const pointField& blockNpoints = blocks[blockNlabel].points();
            const labelList& blockNfaces = blockCells[blockNlabel];

            foundFace = false;
            label blockNfaceLabel;
            for
            (
                blockNfaceLabel = 0;
                blockNfaceLabel < blockNfaces.size();
                blockNfaceLabel++
            )
            {
                if
                (
                    blockFaces[blockNfaces[blockNfaceLabel]]
                 == blockFaces[blockFaceLabel]
                )
                {
                    foundFace = true;
                    break;
                }
            }

            if (!foundFace)
            {
                FatalErrorIn("blockMesh::createMergeList()")
                << "Cannot find merge face for block " << blockNlabel
                    << exit(FatalError);
            };

            const labelListList& blockNfaceFaces =
            blocks[blockNlabel].boundaryPatches()[blockNfaceLabel];

            if (blockPfaceFaces.size() != blockNfaceFaces.size())
            {
                FatalErrorIn("blockMesh::createMergeList()")
                << "Inconsistent number of faces between block pair "
                    << blockPlabel << " and " << blockNlabel
                    << exit(FatalError);
            }


            bool found = false;

            // N-squared point search over all points of all faces of
            // master block over all point of all faces of slave block
            forAll(blockPfaceFaces, blockPfaceFaceLabel)
            {
                const labelList& blockPfaceFacePoints =
                    blockPfaceFaces[blockPfaceFaceLabel];

                labelList& cp = curPairs[blockPfaceFaceLabel];
                cp.setSize(blockPfaceFacePoints.size());
                cp = -1;

                forAll(blockPfaceFacePoints, blockPfaceFacePointLabel)
                {
                    found = false;

                    forAll(blockNfaceFaces, blockNfaceFaceLabel)
                    {
                        const labelList& blockNfaceFacePoints
                        = blockNfaceFaces[blockNfaceFaceLabel];

                        forAll(blockNfaceFacePoints, blockNfaceFacePointLabel)
                        {
                            if
                            (
                                magSqr
                                (
                                    blockPpoints
                                    [
                                        blockPfaceFacePoints
                                        [blockPfaceFacePointLabel]
                                    ]
                                    - blockNpoints
                                    [
                                        blockNfaceFacePoints
                                        [blockNfaceFacePointLabel]
                                    ]
                                ) < sqrMergeTol
                            )
                            {
                                // Found a new pair
                                found = true;

                                cp[blockPfaceFacePointLabel] =
                                blockNfaceFacePoints[blockNfaceFacePointLabel];

                                label PpointLabel =
                                    blockPfaceFacePoints
                                    [
                                        blockPfaceFacePointLabel
                                    ]
                                  + blockOffsets_[blockPlabel];

                                label NpointLabel =
                                    blockNfaceFacePoints
                                    [
                                        blockNfaceFacePointLabel
                                    ]
                                  + blockOffsets_[blockNlabel];

                                label minPN = min(PpointLabel, NpointLabel);

                                if (MergeList[PpointLabel] != -1)
                                {
                                    minPN = min(minPN, MergeList[PpointLabel]);
                                }

                                if (MergeList[NpointLabel] != -1)
                                {
                                    minPN = min(minPN, MergeList[NpointLabel]);
                                }

                                MergeList[PpointLabel] =
                                MergeList[NpointLabel] =
                                minPN;
                            }
                        }
                    }
                }
                forAll(blockPfaceFacePoints, blockPfaceFacePointLabel)
                {
                    if (cp[blockPfaceFacePointLabel] == -1)
                    {
                        FatalErrorIn("blockMesh::createMergeList()")
                            << "Inconsistent point locations between blocks "
                            << blockPlabel << " and " << blockNlabel << nl
                            << "    probably due to inconsistent grading."
                            << exit(FatalError);
                    }
                }
            }
        }
    }


    const faceList::subList blockInternalFaces
    (
        blockFaces,
        topology().nInternalFaces()
    );

    bool changedPointMerge = false;
    label nPasses = 0;

    do
    {
        changedPointMerge = false;
        nPasses++;

        forAll(blockInternalFaces, blockFaceLabel)
        {
            label blockPlabel = faceOwnerBlocks[blockFaceLabel];
            label blockNlabel = faceNeighbourBlocks[blockFaceLabel];

            const labelList& blockPfaces = blockCells[blockPlabel];
            const labelList& blockNfaces = blockCells[blockNlabel];

            const labelListList& curPairs = glueMergePairs[blockFaceLabel];

            bool foundFace = false;
            label blockPfaceLabel;
            for
            (
                blockPfaceLabel = 0;
                blockPfaceLabel < blockPfaces.size();
                blockPfaceLabel++
            )
            {
                if
                (
                    blockFaces[blockPfaces[blockPfaceLabel]]
                 == blockInternalFaces[blockFaceLabel]
                )
                {
                    foundFace = true;
                    break;
                }
            }

            foundFace = false;
            label blockNfaceLabel;
            for
            (
                blockNfaceLabel = 0;
                blockNfaceLabel < blockNfaces.size();
                blockNfaceLabel++
            )
            {
                if
                (
                    blockFaces[blockNfaces[blockNfaceLabel]]
                 == blockInternalFaces[blockFaceLabel]
                )
                {
                    foundFace = true;
                    break;
                }
            }

            const labelListList& blockPfaceFaces =
                blocks[blockPlabel].boundaryPatches()[blockPfaceLabel];

            forAll(blockPfaceFaces, blockPfaceFaceLabel)
            {
                const labelList& blockPfaceFacePoints
                    = blockPfaceFaces[blockPfaceFaceLabel];

                const labelList& cp = curPairs[blockPfaceFaceLabel];

                forAll(blockPfaceFacePoints, blockPfaceFacePointLabel)
                {
                    label PpointLabel =
                        blockPfaceFacePoints[blockPfaceFacePointLabel]
                      + blockOffsets_[blockPlabel];

                    label NpointLabel =
                        cp[blockPfaceFacePointLabel]
                      + blockOffsets_[blockNlabel];

                    if
                    (
                        MergeList[PpointLabel]
                     != MergeList[NpointLabel]
                    )
                    {
                        changedPointMerge = true;

                        MergeList[PpointLabel]
                      = MergeList[NpointLabel]
                      = min
                        (
                            MergeList[PpointLabel],
                            MergeList[NpointLabel]
                        );
                    }
                }
            }
        }
        Info << "." << flush;

        if (nPasses > 100)
        {
            FatalErrorIn("blockMesh::createMergeList()")
                << "Point merging failed after max number of passes."
                << abort(FatalError);
        }
    }
    while (changedPointMerge);
    Info << endl;

    forAll(blockInternalFaces, blockFaceLabel)
    {
        label blockPlabel = faceOwnerBlocks[blockFaceLabel];
        label blockNlabel = faceNeighbourBlocks[blockFaceLabel];

        const labelList& blockPfaces = blockCells[blockPlabel];
        const labelList& blockNfaces = blockCells[blockNlabel];

        const pointField& blockPpoints = blocks[blockPlabel].points();
        const pointField& blockNpoints = blocks[blockNlabel].points();

        bool foundFace = false;
        label blockPfaceLabel;
        for
        (
            blockPfaceLabel = 0;
            blockPfaceLabel < blockPfaces.size();
            blockPfaceLabel++
        )
        {
            if
            (
                blockFaces[blockPfaces[blockPfaceLabel]]
             == blockInternalFaces[blockFaceLabel]
            )
            {
                foundFace = true;
                break;
            }
        }

        if (!foundFace)
        {
            FatalErrorIn("blockMesh::createMergeList()")
                << "Cannot find merge face for block " << blockPlabel
                << exit(FatalError);
        };

        foundFace = false;
        label blockNfaceLabel;
        for
        (
            blockNfaceLabel = 0;
            blockNfaceLabel < blockNfaces.size();
            blockNfaceLabel++
        )
        {
            if
            (
                blockFaces[blockNfaces[blockNfaceLabel]]
             == blockInternalFaces[blockFaceLabel]
            )
            {
                foundFace = true;
                break;
            }
        }

        if (!foundFace)
        {
            FatalErrorIn("blockMesh::createMergeList()")
                << "Cannot find merge face for block " << blockNlabel
                << exit(FatalError);
        };

        const labelListList& blockPfaceFaces =
            blocks[blockPlabel].boundaryPatches()[blockPfaceLabel];

        const labelListList& blockNfaceFaces =
            blocks[blockNlabel].boundaryPatches()[blockNfaceLabel];

        forAll(blockPfaceFaces, blockPfaceFaceLabel)
        {
            const labelList& blockPfaceFacePoints
                = blockPfaceFaces[blockPfaceFaceLabel];

            forAll(blockPfaceFacePoints, blockPfaceFacePointLabel)
            {
                label PpointLabel =
                    blockPfaceFacePoints[blockPfaceFacePointLabel]
                  + blockOffsets_[blockPlabel];

                if (MergeList[PpointLabel] == -1)
                {
                    FatalErrorIn("blockMesh::createMergeList()")
                        << "Unable to merge point "
                        << blockPfaceFacePointLabel
                        << ' ' << blockPpoints[blockPfaceFacePointLabel]
                        << " of face "
                        << blockPfaceLabel
                        << " of block "
                        << blockPlabel
                        << exit(FatalError);
                }
            }
        }

        forAll(blockNfaceFaces, blockNfaceFaceLabel)
        {
            const labelList& blockNfaceFacePoints
                = blockNfaceFaces[blockNfaceFaceLabel];

            forAll(blockNfaceFacePoints, blockNfaceFacePointLabel)
            {
                label NpointLabel =
                    blockNfaceFacePoints[blockNfaceFacePointLabel]
                  + blockOffsets_[blockNlabel];

                if (MergeList[NpointLabel] == -1)
                {
                    FatalErrorIn("blockMesh::createMergeList()")
                        << "unable to merge point "
                        << blockNfaceFacePointLabel
                        << ' ' << blockNpoints[blockNfaceFacePointLabel]
                        << " of face "
                        << blockNfaceLabel
                        << " of block "
                        << blockNlabel
                        << exit(FatalError);
                }
            }
        }
    }


    // sort merge list to return new point label (in new shorter list)
    // given old point label
    label newPointLabel = 0;

    forAll(MergeList, pointLabel)
    {
        if (MergeList[pointLabel] > pointLabel)
        {
            FatalErrorIn("blockMesh::createMergeList()")
                << "ouch" << exit(FatalError);
        }

        if
        (
            (MergeList[pointLabel] == -1)
         || MergeList[pointLabel] == pointLabel
        )
        {
            MergeList[pointLabel] = newPointLabel;
            newPointLabel++;
        }
        else
        {
            MergeList[pointLabel] = MergeList[MergeList[pointLabel]];
        }
    }

    nPoints_ = newPointLabel;


    return MergeList;
}

// ************************************************************************* //
