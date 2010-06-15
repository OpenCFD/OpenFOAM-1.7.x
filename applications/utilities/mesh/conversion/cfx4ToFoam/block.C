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

#include "error.H"

#include "block.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


label block::vtxLabel(label a, label b, label c)
{
    return (a + b*(BlockDef.xDim() + 1)
            + c*(BlockDef.xDim() + 1)*(BlockDef.yDim() + 1));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from description
block::block(const blockDescriptor& definition)
:
    BlockDef(definition),
    Vertices
    (
        ((BlockDef.xDim() + 1)*(BlockDef.yDim() + 1)*(BlockDef.zDim() + 1))
    ),
    Cells
    (
        (BlockDef.xDim()*BlockDef.yDim()*BlockDef.zDim())
    ),
    BoundaryPatches(6)
{
    // create points
    blockPoints();

    // generate internal cells
    blockCells();

    // generate boundary patches
    blockBoundary();
}

// as copy
block::block(const block& original)
:
    BlockDef(original.blockDef()),
    Vertices(original.points()),
    Cells(original.cells()),
    BoundaryPatches(original.boundaryPatches())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

block::~block()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


const blockDescriptor& block::blockDef() const
{
    return BlockDef;
}

const pointField& block::points() const
{
    return Vertices;
}

const labelListList& block::cells() const
{
    return Cells;
}

const labelListListList& block::boundaryPatches() const
{
    return BoundaryPatches;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const block& b)
{
    os << b.Vertices << nl
       << b.Cells << nl
       << b.BoundaryPatches << endl;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

