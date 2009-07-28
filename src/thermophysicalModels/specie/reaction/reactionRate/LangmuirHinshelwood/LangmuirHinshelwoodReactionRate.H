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

Class
    Foam::LangmuirHinshelwoodReactionRate

Description
    Power series reaction rate.

SourceFiles
    LangmuirHinshelwoodReactionRateI.H

\*---------------------------------------------------------------------------*/

#ifndef LangmuirHinshelwoodReactionRate_H
#define LangmuirHinshelwoodReactionRate_H

#include "scalarField.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class LangmuirHinshelwoodReactionRate Declaration
\*---------------------------------------------------------------------------*/

class LangmuirHinshelwoodReactionRate
{
    // Private data

        static const label n_ = 5;
        scalar A_[n_];
        scalar Ta_[n_];

        label co_;
        label c3h6_;
        label no_;


public:

    // Constructors

        //- Construct from components
        inline LangmuirHinshelwoodReactionRate
        (
            const scalar A[],
            const scalar Ta[],
            const label co,
            const label c3h6,
            const label no
        );

        //- Construct from Istream
        inline LangmuirHinshelwoodReactionRate
        (
            const speciesTable& species,
            Istream& is
        );


    // Member Functions

        //- Return the type name
        static word type()
        {
            return "LangmuirHinshelwood";
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
            const LangmuirHinshelwoodReactionRate&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LangmuirHinshelwoodReactionRateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
