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

#include "directMappedWallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directMappedWallPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, directMappedWallPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        directMappedWallPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::directMappedWallPolyPatch::directMappedWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(name, size, start, index, bm),
    directMappedPatchBase(static_cast<const polyPatch&>(*this))
{}


Foam::directMappedWallPolyPatch::directMappedWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const directMappedPatchBase::sampleMode mode,
    const word& samplePatch,
    const vectorField& offset,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(name, size, start, index, bm),
    directMappedPatchBase
    (
        static_cast<const polyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::directMappedWallPolyPatch::directMappedWallPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const directMappedPatchBase::sampleMode mode,
    const word& samplePatch,
    const vector& offset,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(name, size, start, index, bm),
    directMappedPatchBase
    (
        static_cast<const polyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::directMappedWallPolyPatch::directMappedWallPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(name, dict, index, bm),
    directMappedPatchBase(*this, dict)
{}


Foam::directMappedWallPolyPatch::directMappedWallPolyPatch
(
    const directMappedWallPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    wallPolyPatch(pp, bm),
    directMappedPatchBase(*this, pp)
{}


Foam::directMappedWallPolyPatch::directMappedWallPolyPatch
(
    const directMappedWallPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    wallPolyPatch(pp, bm, index, newSize, newStart),
    directMappedPatchBase(*this, pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::directMappedWallPolyPatch::~directMappedWallPolyPatch()
{
    directMappedPatchBase::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Initialise the calculation of the patch geometry
void Foam::directMappedWallPolyPatch::initGeometry()
{
    wallPolyPatch::initGeometry();
    directMappedPatchBase::clearOut();
}

//- Calculate the patch geometry
void Foam::directMappedWallPolyPatch::calcGeometry()
{
    wallPolyPatch::calcGeometry();
    directMappedPatchBase::clearOut();
}

//- Initialise the patches for moving points
void Foam::directMappedWallPolyPatch::initMovePoints(const pointField& p)
{
    wallPolyPatch::initMovePoints(p);
    directMappedPatchBase::clearOut();
}

//- Correct patches after moving points
void Foam::directMappedWallPolyPatch::movePoints(const pointField& p)
{
    wallPolyPatch::movePoints(p);
    directMappedPatchBase::clearOut();
}

//- Initialise the update of the patch topology
void Foam::directMappedWallPolyPatch::initUpdateMesh()
{
    wallPolyPatch::initUpdateMesh();
    directMappedPatchBase::clearOut();
}

//- Update of the patch topology
void Foam::directMappedWallPolyPatch::updateMesh()
{
    wallPolyPatch::updateMesh();
    directMappedPatchBase::clearOut();
}


void Foam::directMappedWallPolyPatch::write(Ostream& os) const
{
    wallPolyPatch::write(os);
    directMappedPatchBase::write(os);
}


// ************************************************************************* //
