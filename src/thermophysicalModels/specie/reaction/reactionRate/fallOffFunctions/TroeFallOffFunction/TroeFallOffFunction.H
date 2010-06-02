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
    Foam::TroeFallOffFunction

Description
    The Troe fall-off function

SourceFiles
    TroeFallOffFunctionI.H

\*---------------------------------------------------------------------------*/

#ifndef TroeFallOffFunction_H
#define TroeFallOffFunction_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class TroeFallOffFunction;
Ostream& operator<<(Ostream&, const TroeFallOffFunction&);


/*---------------------------------------------------------------------------*\
                           Class TroeFallOffFunction Declaration
\*---------------------------------------------------------------------------*/

class TroeFallOffFunction
{
    // Private data

        scalar alpha_;
        scalar Tsss_, Ts_, Tss_;


public:

    // Constructors

        //- Construct from components
        inline TroeFallOffFunction
        (
            const scalar alpha,
            const scalar Tsss,
            const scalar Ts,
            const scalar Tss
        );

        //- Construct from Istream
        inline TroeFallOffFunction(Istream&);


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "Troe";
        }

        inline scalar operator()
        (
            const scalar T,
            const scalar Pr
        ) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const TroeFallOffFunction&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TroeFallOffFunctionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
