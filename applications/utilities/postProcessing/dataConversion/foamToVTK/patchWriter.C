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

#include "patchWriter.H"
#include "writeFuns.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::patchWriter::patchWriter
(
    const vtkMesh& vMesh,
    const bool binary,
    const bool nearCellValue,
    const fileName& fName,
    const labelList& patchIDs
)
:
    vMesh_(vMesh),
    binary_(binary),
    nearCellValue_(nearCellValue),
    fName_(fName),
    patchIDs_(patchIDs),
    os_(fName.c_str())
{
    const fvMesh& mesh = vMesh_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Write header
    if (patchIDs_.size() == 1)
    {
        writeFuns::writeHeader(os_, binary_, patches[patchIDs_[0]].name());
    }
    else
    {
        writeFuns::writeHeader(os_, binary_, "patches");
    }
    os_ << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Write topology
    nPoints_ = 0;
    nFaces_ = 0;
    label nFaceVerts = 0;

    forAll(patchIDs_, i)
    {
        const polyPatch& pp = patches[patchIDs_[i]];

        nPoints_ += pp.nPoints();
        nFaces_ += pp.size();

        forAll(pp, faceI)
        {
            nFaceVerts += pp[faceI].size() + 1;
        }
    }

    os_ << "POINTS " << nPoints_ << " float" << std::endl;

    DynamicList<floatScalar> ptField(3*nPoints_);

    forAll(patchIDs_, i)
    {
        const polyPatch& pp = patches[patchIDs_[i]];

        writeFuns::insert(pp.localPoints(), ptField);
    }
    writeFuns::write(os_, binary_, ptField);

    os_ << "CELLS " << nFaces_ << ' ' << nFaceVerts
        << std::endl;

    DynamicList<label> vertLabels(nFaceVerts);
    DynamicList<label> faceTypes(nFaceVerts);

    label offset = 0;

    forAll(patchIDs_, i)
    {
        const polyPatch& pp = patches[patchIDs_[i]];

        forAll(pp, faceI)
        {
            const face& f = pp.localFaces()[faceI];

            const label fSize = f.size();
            vertLabels.append(fSize);

            writeFuns::insert(f + offset, vertLabels);

            if (fSize == 3)
            {
                faceTypes.append(vtkTopo::VTK_TRIANGLE);
            }
            else if (fSize == 4)
            {
                faceTypes.append(vtkTopo::VTK_QUAD);
            }
            else
            {
                faceTypes.append(vtkTopo::VTK_POLYGON);
            }
        }
        offset += pp.nPoints();
    }
    writeFuns::write(os_, binary_, vertLabels);

    os_ << "CELL_TYPES " << nFaces_ << std::endl;

    writeFuns::write(os_, binary_, faceTypes);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchWriter::writePatchIDs()
{
    const fvMesh& mesh = vMesh_.mesh();

    DynamicList<floatScalar> fField(nFaces_);

    os_ << "patchID 1 " << nFaces_ << " float" << std::endl;

    forAll(patchIDs_, i)
    {
        label patchI = patchIDs_[i];

        const polyPatch& pp = mesh.boundaryMesh()[patchI];

        if (!isA<emptyPolyPatch>(pp))
        {
            writeFuns::insert(scalarField(pp.size(), patchI), fField);
        }
    }
    writeFuns::write(os_, binary_, fField);
}


// ************************************************************************* //
