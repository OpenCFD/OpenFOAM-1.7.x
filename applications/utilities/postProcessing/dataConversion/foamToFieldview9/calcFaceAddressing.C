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

#include "calcFaceAddressing.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// Returns the face labels of the shape in an order consistent with the
// shape.
labelList calcFaceAddressing
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
                //Info<< "match:" << modelFace
                //    << "  to " << vertLabels << endl;
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


// ************************************************************************* //
