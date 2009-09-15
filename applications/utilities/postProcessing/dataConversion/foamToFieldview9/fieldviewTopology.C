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

\*---------------------------------------------------------------------------*/

#include "fieldviewTopology.H"
#include "polyMesh.H"
#include "cellShape.H"
#include "cellModeller.H"
#include "wallPolyPatch.H"
#include "symmetryPolyPatch.H"


#include "fv_reader_tags.h"

extern "C"
{
    unsigned int fv_encode_elem_header(int elem_type, int wall_info[]);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::fieldviewTopology::calcFaceAddressing
(
    const faceList& allFaces,   // faces given faceLabels
    const cellShape& shape,
    const labelList& faces,     // faceLabels for given cell
    const label cellI
)
{
    // return value. 
    labelList shapeToMesh(shape.nFaces(), -1);

    const faceList modelFaces(shape.faces());

    // Loop over all faces of cellShape
    forAll(modelFaces, cellFaceI)
    {
        // face (vertex list)
        const face& modelFace = modelFaces[cellFaceI];

        // Loop over all face labels
        forAll(faces, faceI)
        {
            const face& vertLabels = allFaces[faces[faceI]];

            if (vertLabels == modelFace)
            {
                shapeToMesh[cellFaceI] = faces[faceI];
                break;
            }
        }

        if (shapeToMesh[cellFaceI] == -1)
        {
            FatalErrorIn("foamToFieldview : calcFaceAddressing")
                << "calcFaceAddressing : can't match face to shape.\n"
                << "    shape face:" << modelFace << endl
                << "    face labels:" << faces << endl
                << "    cellI:" << cellI << endl;

            FatalError << "Faces consist of vertices:" << endl;
            forAll(faces, faceI)
            {
                FatalError
                    << "    face:" << faces[faceI]
                    << allFaces[faces[faceI]] << endl;
            }
            FatalError << exit(FatalError);
        }
    }
    return shapeToMesh;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::fieldviewTopology::fieldviewTopology
(
    const polyMesh& mesh,
    const bool setWallInfo
)
:
    hexLabels_((1+8)*mesh.nCells()),
    prismLabels_((1+6)*mesh.nCells()),
    pyrLabels_((1+5)*mesh.nCells()),
    tetLabels_((1+4)*mesh.nCells()),
    nPoly_(0),
    quadFaceLabels_(mesh.boundaryMesh().size()),
    nPolyFaces_(mesh.boundaryMesh().size())
{
    // Mark all faces that are to be seen as wall for particle
    // tracking and all cells that use one or more of these walls

    labelList wallFace(mesh.nFaces(), NOT_A_WALL);
    boolList wallCell(mesh.nCells(), false);

    if (setWallInfo)
    {
        forAll (mesh.boundaryMesh(), patchI)
        {
            const polyPatch& currPatch = mesh.boundaryMesh()[patchI];
            if 
            (
                isA<wallPolyPatch>(currPatch)
             || isA<symmetryPolyPatch>(currPatch)
            )
            {
                forAll(currPatch, patchFaceI)
                {
                    label meshFaceI = currPatch.start() + patchFaceI;

                    wallFace[meshFaceI] = A_WALL;
                    wallCell[mesh.faceOwner()[meshFaceI]] = true;
                }
            }
        }
    }



    const cellModel& tet = *(cellModeller::lookup("tet"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& wedge = *(cellModeller::lookup("wedge"));
    const cellModel& tetWedge = *(cellModeller::lookup("tetWedge"));
    const cellModel& hex = *(cellModeller::lookup("hex"));

    // Pre calculate headers for cells not on walls
    labelList notWallFlags(6, NOT_A_WALL);
    label tetNotWall = fv_encode_elem_header
    (
        FV_TET_ELEM_ID, notWallFlags.begin()
    );
    label pyrNotWall = fv_encode_elem_header
    (
        FV_PYRA_ELEM_ID, notWallFlags.begin()
    );
    label prismNotWall = fv_encode_elem_header
    (
        FV_PRISM_ELEM_ID, notWallFlags.begin()
    );
    label hexNotWall = fv_encode_elem_header
    (
        FV_HEX_ELEM_ID, notWallFlags.begin()
    );

    // Some aliases
    const cellList& cellFaces = mesh.cells();
    const cellShapeList& cellShapes = mesh.cellShapes();


    label hexi = 0;
    label prismi = 0;
    label pyri = 0;
    label teti = 0;

    const faceList& allFaces = mesh.faces();

    labelList wallFlags(6);
    forAll(cellShapes, celli)
    {
        const cellShape& cellShape = cellShapes[celli];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == tet)
        {
            if (!wallCell[celli])
            {
                tetLabels_[teti++] = tetNotWall;
            }
            else
            {
                labelList modelToMesh = calcFaceAddressing
                (
                    allFaces, cellShape, cellFaces[celli], celli
                );

                wallFlags[0] = wallFace[modelToMesh[0]];
                wallFlags[1] = wallFace[modelToMesh[1]];
                wallFlags[2] = wallFace[modelToMesh[2]];
                wallFlags[3] = wallFace[modelToMesh[3]];

                tetLabels_[teti++] = fv_encode_elem_header
                (
                    FV_TET_ELEM_ID, wallFlags.begin()
                );
            }

            tetLabels_[teti++] = cellShape[0] + 1;
            tetLabels_[teti++] = cellShape[1] + 1;
            tetLabels_[teti++] = cellShape[2] + 1;
            tetLabels_[teti++] = cellShape[3] + 1;
        }
        else if (cellModel == pyr)
        {
            if (!wallCell[celli])
            {
                pyrLabels_[pyri++] = pyrNotWall;
            }
            else
            {
                labelList modelToMesh = calcFaceAddressing
                (
                    allFaces, cellShape, cellFaces[celli], celli
                );

                wallFlags[0] = wallFace[modelToMesh[0]];
                wallFlags[1] = wallFace[modelToMesh[3]];
                wallFlags[2] = wallFace[modelToMesh[2]];
                wallFlags[3] = wallFace[modelToMesh[1]];
                wallFlags[4] = wallFace[modelToMesh[4]];

                pyrLabels_[pyri++] = fv_encode_elem_header
                (
                    FV_PYRA_ELEM_ID, wallFlags.begin()
                );
            }

            pyrLabels_[pyri++] = cellShape[0] + 1;
            pyrLabels_[pyri++] = cellShape[1] + 1;
            pyrLabels_[pyri++] = cellShape[2] + 1;
            pyrLabels_[pyri++] = cellShape[3] + 1;
            pyrLabels_[pyri++] = cellShape[4] + 1;
        }
        else if (cellModel == prism)
        {
            if (!wallCell[celli])
            {
                prismLabels_[prismi++] = prismNotWall;
            }
            else
            {
                labelList modelToMesh = calcFaceAddressing
                (
                    allFaces, cellShape, cellFaces[celli], celli
                );

                wallFlags[0] = wallFace[modelToMesh[4]];
                wallFlags[1] = wallFace[modelToMesh[2]];
                wallFlags[2] = wallFace[modelToMesh[3]];
                wallFlags[3] = wallFace[modelToMesh[0]];
                wallFlags[4] = wallFace[modelToMesh[1]];

                prismLabels_[prismi++] = fv_encode_elem_header
                (
                    FV_PRISM_ELEM_ID, wallFlags.begin()
                );
            }

            prismLabels_[prismi++] = cellShape[0] + 1;
            prismLabels_[prismi++] = cellShape[3] + 1;
            prismLabels_[prismi++] = cellShape[4] + 1;
            prismLabels_[prismi++] = cellShape[1] + 1;
            prismLabels_[prismi++] = cellShape[5] + 1;
            prismLabels_[prismi++] = cellShape[2] + 1;
        }
        else if (cellModel == tetWedge)
        {
            // Treat as prism with collapsed edge
            if (!wallCell[celli])
            {
                prismLabels_[prismi++] = prismNotWall;
            }
            else
            {
                labelList modelToMesh = calcFaceAddressing
                (
                    allFaces, cellShape, cellFaces[celli], celli
                );

                wallFlags[0] = wallFace[modelToMesh[1]];
                wallFlags[1] = wallFace[modelToMesh[2]];
                wallFlags[2] = wallFace[modelToMesh[3]];
                wallFlags[3] = wallFace[modelToMesh[0]];
                wallFlags[4] = wallFace[modelToMesh[3]];

                prismLabels_[prismi++] = fv_encode_elem_header
                (
                    FV_PRISM_ELEM_ID, wallFlags.begin()
                );
            }

            prismLabels_[prismi++] = cellShape[0] + 1;
            prismLabels_[prismi++] = cellShape[3] + 1;
            prismLabels_[prismi++] = cellShape[4] + 1;
            prismLabels_[prismi++] = cellShape[1] + 1;
            prismLabels_[prismi++] = cellShape[4] + 1;
            prismLabels_[prismi++] = cellShape[2] + 1;
        }
        else if (cellModel == wedge)
        {
            if (!wallCell[celli])
            {
                hexLabels_[hexi++] = hexNotWall;
            }
            else
            {
                labelList modelToMesh = calcFaceAddressing
                (
                    allFaces, cellShape, cellFaces[celli], celli
                );

                wallFlags[0] = wallFace[modelToMesh[2]];
                wallFlags[1] = wallFace[modelToMesh[3]];
                wallFlags[2] = wallFace[modelToMesh[0]];
                wallFlags[3] = wallFace[modelToMesh[1]];
                wallFlags[4] = wallFace[modelToMesh[4]];
                wallFlags[5] = wallFace[modelToMesh[5]];

                hexLabels_[hexi++] = fv_encode_elem_header
                (
                    FV_HEX_ELEM_ID, wallFlags.begin()
                );
            }
            hexLabels_[hexi++] = cellShape[0] + 1;
            hexLabels_[hexi++] = cellShape[1] + 1;
            hexLabels_[hexi++] = cellShape[0] + 1;
            hexLabels_[hexi++] = cellShape[2] + 1;
            hexLabels_[hexi++] = cellShape[3] + 1;
            hexLabels_[hexi++] = cellShape[4] + 1;
            hexLabels_[hexi++] = cellShape[6] + 1;
            hexLabels_[hexi++] = cellShape[5] + 1;
        }
        else if (cellModel == hex)
        {
            if (!wallCell[celli])
            {
                hexLabels_[hexi++] = hexNotWall;
            }
            else
            {
                labelList modelToMesh = calcFaceAddressing
                (
                    allFaces, cellShape, cellFaces[celli], celli
                );

                wallFlags[0] = wallFace[modelToMesh[0]];
                wallFlags[1] = wallFace[modelToMesh[1]];
                wallFlags[2] = wallFace[modelToMesh[4]];
                wallFlags[3] = wallFace[modelToMesh[5]];
                wallFlags[4] = wallFace[modelToMesh[2]];
                wallFlags[5] = wallFace[modelToMesh[3]];

                hexLabels_[hexi++] = fv_encode_elem_header
                (
                    FV_HEX_ELEM_ID, wallFlags.begin()
                );
            }
            hexLabels_[hexi++] = cellShape[0] + 1;
            hexLabels_[hexi++] = cellShape[1] + 1;
            hexLabels_[hexi++] = cellShape[3] + 1;
            hexLabels_[hexi++] = cellShape[2] + 1;
            hexLabels_[hexi++] = cellShape[4] + 1;
            hexLabels_[hexi++] = cellShape[5] + 1;
            hexLabels_[hexi++] = cellShape[7] + 1;
            hexLabels_[hexi++] = cellShape[6] + 1;
        }
        else
        {
            nPoly_++;
        }
    }

    hexLabels_.setSize(hexi);
    prismLabels_.setSize(prismi);
    pyrLabels_.setSize(pyri);
    tetLabels_.setSize(teti);


    //
    // Patches
    //
    forAll(mesh.boundaryMesh(), patchI)
    {
        const polyPatch& patchFaces = mesh.boundaryMesh()[patchI];

        labelList& faceLabels = quadFaceLabels_[patchI];

        // Faces, each 4 labels. Size big enough
        faceLabels.setSize(patchFaces.size()*4);

        label labelI = 0;

        forAll(patchFaces, faceI)
        {
            const face& patchFace = patchFaces[faceI];

            if (patchFace.size() == 3)
            {
                faceLabels[labelI++] = patchFace[0] + 1;
                faceLabels[labelI++] = patchFace[1] + 1;
                faceLabels[labelI++] = patchFace[2] + 1;
                faceLabels[labelI++] = 0;   // Fieldview:triangle definition
            }
            else if (patchFace.size() == 4)
            {
                faceLabels[labelI++] = patchFace[0] + 1;
                faceLabels[labelI++] = patchFace[1] + 1;
                faceLabels[labelI++] = patchFace[2] + 1;
                faceLabels[labelI++] = patchFace[3] + 1;
            } 
        }

        faceLabels.setSize(labelI);

        label nFaces = labelI/4;

        nPolyFaces_[patchI] = patchFaces.size() - nFaces;
    }
}


// ************************************************************************* //
