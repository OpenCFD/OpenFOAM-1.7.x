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


\*---------------------------------------------------------------------------*/

#include "blockMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::faceList Foam::blockMesh::createPatchFaces
(
    const polyPatch& patchTopologyFaces
)
{
    blockMesh& blocks = *this;

    labelList blockLabels = patchTopologyFaces.polyPatch::faceCells();

    label nFaces=0;

    forAll(patchTopologyFaces, patchTopologyFaceLabel)
    {
        label blockLabel = blockLabels[patchTopologyFaceLabel];

        faceList blockFaces
        (
            blocks[blockLabel].blockDef().blockShape().faces()
        );

        forAll(blockFaces, blockFaceLabel)
        {
            if
            (
                blockFaces[blockFaceLabel]
             == patchTopologyFaces[patchTopologyFaceLabel]
            )
            {
                nFaces +=
                    blocks[blockLabel].boundaryPatches()[blockFaceLabel].size();
            }
        }
    }


    faceList patchFaces(nFaces);
    face quadFace(4);
    label faceLabel = 0;

    forAll(patchTopologyFaces, patchTopologyFaceLabel)
    {
        label blockLabel = blockLabels[patchTopologyFaceLabel];

        faceList blockFaces
        (
            blocks[blockLabel].blockDef().blockShape().faces()
        );

        forAll(blockFaces, blockFaceLabel)
        {
            if
            (
                blockFaces[blockFaceLabel]
                == patchTopologyFaces[patchTopologyFaceLabel]
            )
            {
                const labelListList& blockPatchFaces =
                    blocks[blockLabel].boundaryPatches()[blockFaceLabel];

                forAll(blockPatchFaces, blockFaceLabel)
                {
                    // Lookup the face points
                    // and collapse duplicate point labels

                    quadFace[0] =
                        mergeList_
                        [
                            blockPatchFaces[blockFaceLabel][0]
                          + blockOffsets_[blockLabel]
                        ];

                    label nUnique = 1;

                    for
                    (
                        label facePointLabel = 1;
                        facePointLabel < 4;
                        facePointLabel++
                    )
                    {
                        quadFace[nUnique] =
                            mergeList_
                            [
                                blockPatchFaces[blockFaceLabel][facePointLabel]
                              + blockOffsets_[blockLabel]
                            ];

                        if (quadFace[nUnique] != quadFace[nUnique-1])
                        {
                            nUnique++;
                        }
                    }

                    if (quadFace[nUnique-1] == quadFace[0])
                    {
                        nUnique--;
                    }

                    if (nUnique == 4)
                    {
                        patchFaces[faceLabel++] = quadFace;
                    }
                    else if (nUnique == 3)
                    {
                        patchFaces[faceLabel++] = face
                        (
                            labelList::subList(quadFace, 3)
                        );
                    }
                    // else the face has collapsed to an edge or point
                }
            }
        }
    }

    patchFaces.setSize(faceLabel);

    return patchFaces;
}


Foam::faceListList Foam::blockMesh::createPatches()
{
    Info<< "\nCreating patches\n";

    const polyPatchList& patchTopologies = topology().boundaryMesh();
    faceListList patches(patchTopologies.size());

    forAll(patchTopologies, patchLabel)
    {
        patches[patchLabel] =
            createPatchFaces(patchTopologies[patchLabel]);
    }

    return patches;
}

// ************************************************************************* //
