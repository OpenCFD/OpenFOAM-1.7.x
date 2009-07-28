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

#include "CarnahanStarlingRadial.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CarnahanStarlingRadial, 0);

    addToRunTimeSelectionTable
    (
        radialModel,
        CarnahanStarlingRadial,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CarnahanStarlingRadial::CarnahanStarlingRadial(const dictionary& dict)
:
    radialModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::CarnahanStarlingRadial::~CarnahanStarlingRadial()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::CarnahanStarlingRadial::g0
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax
) const
{

    return 
        1.0/(1.0 - alpha)
      + 3.0*alpha/(2.0*sqr(1.0 - alpha))
      + sqr(alpha)/(2.0*pow(1.0 - alpha, 3));
}


Foam::tmp<Foam::volScalarField> Foam::CarnahanStarlingRadial::g0prime
(
    const volScalarField& alpha,
    const dimensionedScalar& alphaMax
) const
{
    return 
        - alpha/sqr(1.0 - alpha) 
        + (3.0*(1.0 - alpha) + 6.0*sqr(alpha))/(2.0*(1.0 - alpha))
        + (2.0*alpha*(1.0 - alpha) + 3.0*pow(alpha, 3))
         /(2.0*pow(1.0 - alpha, 4));
}


// ************************************************************************* //
