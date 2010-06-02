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

Class
    Foam::thirdBodyArrheniusReactionRate

Description
    Arrhenius reaction rate enhanced by third-body interation.

SourceFiles
    thirdBodyArrheniusReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef thirdBodyArrheniusReactionRate_H
#define thirdBodyArrheniusReactionRate_H

#include "ArrheniusReactionRate.H"
#include "thirdBodyEfficiencies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class thirdBodyArrheniusReactionRate Declaration
\*---------------------------------------------------------------------------*/

class thirdBodyArrheniusReactionRate
:
    ArrheniusReactionRate
{
    // Private data

        thirdBodyEfficiencies thirdBodyEfficiencies_;


public:

    // Constructors

        //- Construct from components
        inline thirdBodyArrheniusReactionRate
        (
            const scalar A,
            const scalar beta,
            const scalar Ta,
            const thirdBodyEfficiencies& tbes
        );

        //- Construct from Istream
        inline thirdBodyArrheniusReactionRate
        (
            const speciesTable& species,
            Istream& is
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "thirdBodyArrhenius";
        }

        inline scalar operator()
        (
            const scalar T,
            const scalar p,
            const scalarField& c
        ) const;


    // Ostream Operator

        inline friend Ostream& operator<<
        (
            Ostream&,
            const thirdBodyArrheniusReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thirdBodyArrheniusReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
