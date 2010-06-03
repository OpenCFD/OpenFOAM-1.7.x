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

#include "solid.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solid, 0);
    defineRunTimeSelectionTable(solid,);
    defineRunTimeSelectionTable(solid, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solid::solid
(
    scalar rho,
    scalar cp,
    scalar K,
    scalar Hf,
    scalar emissivity
)
:
    rho_(rho),
    cp_(cp),
    K_(K),
    Hf_(Hf),
    emissivity_(emissivity)
{}


Foam::solid::solid(Istream& is)
:
    rho_(readScalar(is)),
    cp_(readScalar(is)),
    K_(readScalar(is)),
    Hf_(readScalar(is)),
    emissivity_(readScalar(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solid::writeData(Ostream& os) const
{
    os  << rho_ << token::SPACE
        << cp_ << token::SPACE
        << K_ << token::SPACE
        << Hf_ << token::SPACE
        << emissivity_;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const solid& s)
{
    s.writeData(os);
    return os;
}


// ************************************************************************* //
