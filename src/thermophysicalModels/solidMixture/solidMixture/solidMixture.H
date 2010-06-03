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

Class
    Foam::solidMixture

Description
    Foam::solidMixture

\*---------------------------------------------------------------------------*/

#ifndef solidMixture_H
#define solidMixture_H

#include "scalarField.H"
#include "PtrList.H"
#include "solid.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class solidMixture Declaration
\*---------------------------------------------------------------------------*/

class solidMixture
{

    // Private data

        //- The names of the solids
        List<word> components_;

        //- The solid properties
        PtrList<solid> properties_;


public:

    // Constructors

        //- Construct from dictionary
        solidMixture
        (
            const dictionary&
        );


    // Selectors

        //- Select construct from dictionary
        static autoPtr<solidMixture> New
        (
            const dictionary&
        );


    // Member Functions

        //- Return the sold names
        inline const List<word>& components() const
        {
            return components_;
        }

        //- Return the solid properties
        inline const PtrList<solid>& properties() const
        {
            return properties_;
        }

        //- Returns the mass fractions, given mole fractions
        scalarField Y(const scalarField& X) const;

        //- Returns the mole fractions, given mass fractions
        scalarField X(const scalarField& Y) const;

        //- Calculate the mixture density [kg/m^3] as a function of
        //  volume fractions
        scalar rho
        (
            const scalarField& X
        ) const;

        //- Calculate the mixture heat capacity [J/(kg K)] as a function
        //  of mass fractions
        scalar cp
        (
            const scalarField& Y
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
