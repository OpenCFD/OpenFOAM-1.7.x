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

#include "perfectGas.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Construct from components
inline perfectGas::perfectGas
(
    const specie& sp
)
:
    specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct as named copy
inline perfectGas::perfectGas(const word& name, const perfectGas& pg)
:
    specie(name, pg)
{}


// Construct and return a clone
inline autoPtr<perfectGas> perfectGas::clone() const
{
    return autoPtr<perfectGas>(new perfectGas(*this));
}


// Selector from Istream
inline autoPtr<perfectGas> perfectGas::New(Istream& is)
{
    return autoPtr<perfectGas>(new perfectGas(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return density [kg/m^3]
inline scalar perfectGas::rho(scalar p, scalar T) const
{
    return p/(R()*T);
}

//- Return compressibility rho/p [s^2/m^2]
inline scalar perfectGas::psi(scalar, scalar T) const
{
    return 1.0/(R()*T);
}

//- Return compression factor []
inline scalar perfectGas::Z(scalar, scalar) const
{
    return 1.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void perfectGas::operator+=(const perfectGas& pg)
{
    specie::operator+=(pg);
}

inline void perfectGas::operator-=(const perfectGas& pg)
{
    specie::operator-=(pg);
}

inline void perfectGas::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline perfectGas operator+
(
    const perfectGas& pg1,
    const perfectGas& pg2
)
{
    return perfectGas
    (
        static_cast<const specie&>(pg1)
      + static_cast<const specie&>(pg2)
    );
}


inline perfectGas operator-
(
    const perfectGas& pg1,
    const perfectGas& pg2
)
{
    return perfectGas
    (
        static_cast<const specie&>(pg1)
      - static_cast<const specie&>(pg2)
    );
}


inline perfectGas operator*
(
    const scalar s,
    const perfectGas& pg
)
{
    return perfectGas(s*static_cast<const specie&>(pg));
}


inline perfectGas operator==
(
    const perfectGas& pg1,
    const perfectGas& pg2
)
{
    return pg2 - pg1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
