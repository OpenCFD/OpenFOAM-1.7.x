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

#include "patchIdentifier.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchIdentifier::patchIdentifier
(
    const word& name,
    const label index,
    const word& physicalType
)
:
    name_(name),
    index_(index),
    physicalType_(physicalType)
{}


Foam::patchIdentifier::patchIdentifier
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    name_(name),
    index_(index)
{
    dict.readIfPresent("physicalType", physicalType_);
}


Foam::patchIdentifier::patchIdentifier
(
    const patchIdentifier& p,
    const label index
)
:
    name_(p.name_),
    index_(index),
    physicalType_(p.physicalType_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchIdentifier::~patchIdentifier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchIdentifier::write(Ostream& os) const
{
    if (physicalType_.size())
    {
        os.writeKeyword("physicalType") << physicalType_
            << token::END_STATEMENT << nl;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const patchIdentifier& pi)
{
    pi.write(os);
    os.check("Ostream& operator<<(Ostream&, const patchIdentifier&)");
    return os;
}


// ************************************************************************* //
