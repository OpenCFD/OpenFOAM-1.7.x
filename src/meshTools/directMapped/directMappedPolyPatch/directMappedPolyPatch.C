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

#include "directMappedPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directMappedPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, directMappedPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, directMappedPolyPatch, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm),
    directMappedPatchBase(static_cast<const polyPatch&>(*this))
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
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
    polyPatch(name, size, start, index, bm),
    directMappedPatchBase
    (
        static_cast<const polyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
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
    polyPatch(name, size, start, index, bm),
    directMappedPatchBase
    (
        static_cast<const polyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, dict, index, bm),
    directMappedPatchBase(*this, dict)
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const directMappedPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    directMappedPatchBase(*this, pp)
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const directMappedPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    directMappedPatchBase(*this, pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::directMappedPolyPatch::~directMappedPolyPatch()
{
    directMappedPatchBase::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Initialise the calculation of the patch geometry
void Foam::directMappedPolyPatch::initGeometry()
{
    polyPatch::initGeometry();
    directMappedPatchBase::clearOut();
}

//- Calculate the patch geometry
void Foam::directMappedPolyPatch::calcGeometry()
{
    polyPatch::calcGeometry();
    directMappedPatchBase::clearOut();
}

//- Initialise the patches for moving points
void Foam::directMappedPolyPatch::initMovePoints(const pointField& p)
{
    polyPatch::initMovePoints(p);
    directMappedPatchBase::clearOut();
}

//- Correct patches after moving points
void Foam::directMappedPolyPatch::movePoints(const pointField& p)
{
    polyPatch::movePoints(p);
    directMappedPatchBase::clearOut();
}

//- Initialise the update of the patch topology
void Foam::directMappedPolyPatch::initUpdateMesh()
{
    polyPatch::initUpdateMesh();
    directMappedPatchBase::clearOut();
}

//- Update of the patch topology
void Foam::directMappedPolyPatch::updateMesh()
{
    polyPatch::updateMesh();
    directMappedPatchBase::clearOut();
}


void Foam::directMappedPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    directMappedPatchBase::write(os);
}


// ************************************************************************* //
