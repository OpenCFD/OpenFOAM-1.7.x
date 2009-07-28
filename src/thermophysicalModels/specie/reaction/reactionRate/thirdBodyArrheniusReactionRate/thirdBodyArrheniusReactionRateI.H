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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
inline thirdBodyArrheniusReactionRate::thirdBodyArrheniusReactionRate
(
    const scalar A,
    const scalar beta,
    const scalar Ta,
    const thirdBodyEfficiencies& tbes
)
:
    ArrheniusReactionRate(A, beta, Ta),
    thirdBodyEfficiencies_(tbes)
{}


//- Construct from Istream
inline thirdBodyArrheniusReactionRate::thirdBodyArrheniusReactionRate
(
    const speciesTable& species,
    Istream& is
)
:
    ArrheniusReactionRate
    (
        species, 
        is.readBegin("thirdBodyArrheniusReactionRate(Istream&)")
    ),
    thirdBodyEfficiencies_(species, is)
{
    is.readEnd("thirdBodyArrheniusReactionRate(Istream&)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar thirdBodyArrheniusReactionRate::operator()
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return 
        thirdBodyEfficiencies_.M(c)
       *ArrheniusReactionRate::operator()(T, p, c);
}


inline Ostream& operator<<(Ostream& os, const thirdBodyArrheniusReactionRate& arr)
{
    os  << token::BEGIN_LIST
        << static_cast<const ArrheniusReactionRate&>(arr)
        << token::SPACE << arr.thirdBodyEfficiencies_
        << token::END_LIST;
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
