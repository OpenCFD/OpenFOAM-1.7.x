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

#include "pointBoundaryMesh.H"
#include "polyBoundaryMesh.H"
#include "facePointPatch.H"
#include "globalPointPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointBoundaryMesh::pointBoundaryMesh
(
    const pointMesh& m,
    const polyBoundaryMesh& basicBdry
)
:
    pointPatchList(basicBdry.size()),
    mesh_(m)
{
    // Set boundary patches
    pointPatchList& Patches = *this;

    forAll(Patches, patchI)
    {
        Patches.set
        (
            patchI,
            facePointPatch::New(basicBdry[patchI], *this).ptr()
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointBoundaryMesh::calcGeometry()
{
    forAll(*this, patchi)
    {
        operator[](patchi).initGeometry();
    }

    forAll(*this, patchi)
    {
        operator[](patchi).calcGeometry();
    }
}


const Foam::globalPointPatch&
Foam::pointBoundaryMesh::globalPatch() const
{
    const pointPatchList& patches = *this;

    forAll (patches, patchI)
    {
        if (isType<globalPointPatch>(patches[patchI]))
        {
            return refCast<const globalPointPatch>(patches[patchI]);
        }
    }

    FatalErrorIn
    (
        "const pointBoundaryMesh::"
        "globalPointPatch& globalPatch() const"
    )   << "patch not found."
        << abort(FatalError);

    // Dummy return
    return refCast<const globalPointPatch>(patches[0]);
}


void Foam::pointBoundaryMesh::movePoints(const pointField& p)
{
    pointPatchList& patches = *this;

    forAll(patches, patchi)
    {
        patches[patchi].initMovePoints(p);
    }

    forAll(patches, patchi)
    {
        patches[patchi].movePoints(p);
    }
}


void Foam::pointBoundaryMesh::updateMesh()
{
    pointPatchList& patches = *this;

    forAll(patches, patchi)
    {
        patches[patchi].initUpdateMesh();
    }

    forAll(patches, patchi)
    {
        patches[patchi].updateMesh();
    }
}


// ************************************************************************* //
