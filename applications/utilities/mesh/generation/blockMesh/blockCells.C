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

Description
    private member of block. Creates cells for the block.

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "block.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::block::blockCells()
{
    label ni = blockDef_.n().x();
    label nj = blockDef_.n().y();
    label nk = blockDef_.n().z();

    label cellNo = 0;

    for (label k = 0; k <= nk - 1; k++)
    {
        for (label j = 0; j <= nj - 1; j++)
        {
            for (label i = 0; i <= ni - 1; i++)
            {
                cells_[cellNo].setSize(8);

                cells_[cellNo][0] =  vtxLabel(i, j, k);
                cells_[cellNo][1] =  vtxLabel(i+1, j, k);
                cells_[cellNo][2] =  vtxLabel(i+1, j+1, k);
                cells_[cellNo][3] =  vtxLabel(i, j+1, k);
                cells_[cellNo][4] =  vtxLabel(i, j, k+1);
                cells_[cellNo][5] =  vtxLabel(i+1, j, k+1);
                cells_[cellNo][6] =  vtxLabel(i+1, j+1, k+1);
                cells_[cellNo][7] =  vtxLabel(i, j+1, k+1);
                cellNo++;
            }
        }
    }
}

// ************************************************************************* //
