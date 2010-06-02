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

#include "mapDistribute.H"
#include "extendedCentredFaceToCellStencil.H"
#include "faceToCellStencil.H"

// Only for access to calcDistributeMap <- needs to be moved out
#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedCentredFaceToCellStencil::extendedCentredFaceToCellStencil
(
    const faceToCellStencil& stencil
)
:
    extendedFaceToCellStencil(stencil.mesh())
{
    stencil_ = stencil;

    // Calculate distribute map (also renumbers elements in stencil)
    mapPtr_ = extendedCellToFaceStencil::calcDistributeMap
    (
        stencil.mesh(),
        stencil.globalNumbering(),
        stencil_
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Per face which elements of the stencil to keep.
void Foam::extendedCentredFaceToCellStencil::compact()
{
    boolList isInStencil(map().constructSize(), false);

    forAll(stencil_, faceI)
    {
        const labelList& stencilCells = stencil_[faceI];

        forAll(stencilCells, i)
        {
            isInStencil[stencilCells[i]] = true;
        }
    }

    mapPtr_().compact(isInStencil);
}


// ************************************************************************* //
