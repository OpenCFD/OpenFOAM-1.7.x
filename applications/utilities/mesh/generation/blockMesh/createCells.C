/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "error.H"

#include "blockMesh.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::cellShapeList Foam::blockMesh::createCells()
{
    Info<< nl << "Creating cells" << endl;

    PtrList<cellShape> cells(nCells_);

    blockMesh& blocks = *this;

    const cellModel& hex = *(cellModeller::lookup("hex"));

    label cellLabel = 0;

    forAll(blocks, blockLabel)
    {
        const labelListList& blockCells = blocks[blockLabel].cells();

        forAll(blockCells, blockCellLabel)
        {
            labelList cellPoints(blockCells[blockCellLabel].size());

            forAll(cellPoints, cellPointLabel)
            {
                cellPoints[cellPointLabel] =
                    mergeList_
                    [
                        blockCells[blockCellLabel][cellPointLabel]
                      + blockOffsets_[blockLabel]
                    ];
            }

            // Construct collapsed cell and all to list
            cells.set(cellLabel, new cellShape(hex, cellPoints, true));

            cellLabel++;
        }
    }

    return cells;
}

// ************************************************************************* //
