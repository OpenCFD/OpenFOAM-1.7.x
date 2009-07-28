/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "searchableSurface.H"
#include "long.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchableSurface, 0);
defineRunTimeSelectionTable(searchableSurface, dict);


// Construct named object from dictionary
autoPtr<searchableSurface> searchableSurface::New
(
    const word& searchableSurfaceType,
    const IOobject& io,
    const dictionary& dict
)
{
    dictConstructorTable::iterator cstrIter =
        dictConstructorTablePtr_
            ->find(searchableSurfaceType);

    if (cstrIter == dictConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "searchableSurface::New(const word&, const word&"
            ", const IOobject&, const dictionary&)"
        )   << "Unknown searchableSurface type " << searchableSurfaceType
            << endl << endl
            << "Valid searchableSurface types : " << endl
            << dictConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<searchableSurface>(cstrIter()(io, dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurface::searchableSurface(const IOobject& io)
:
    regIOobject(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurface::~searchableSurface()
{}


} // End namespace Foam


// ************************************************************************* //
