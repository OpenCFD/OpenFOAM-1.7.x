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

#include "boundaryToFace.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(boundaryToFace, 0);

addToRunTimeSelectionTable(topoSetSource, boundaryToFace, word);

addToRunTimeSelectionTable(topoSetSource, boundaryToFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::boundaryToFace::usage_
(
    boundaryToFace::typeName,
    "\n    Usage: boundaryToFace\n\n"
    "    Select all boundary faces\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boundaryToFace::combine(topoSet& set, const bool add) const
{
    for
    (
        label faceI = mesh().nInternalFaces();
        faceI < mesh().nFaces();
        faceI++
    )
    {
        addOrDelete(set, faceI, add);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::boundaryToFace::boundaryToFace(const polyMesh& mesh)
:
    topoSetSource(mesh)
{}


// Construct from dictionary
Foam::boundaryToFace::boundaryToFace(const polyMesh& mesh, const dictionary&)
:
    topoSetSource(mesh)
{}


// Construct from Istream
Foam::boundaryToFace::boundaryToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryToFace::~boundaryToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundaryToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all boundary faces ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all boundary faces ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
