/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    Foam::energyScalingFunctions::shiftedForce

Description


SourceFiles
    shiftedForce.C

\*---------------------------------------------------------------------------*/

#ifndef shiftedForce_H
#define shiftedForce_H

#include "energyScalingFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace energyScalingFunctions
{

/*---------------------------------------------------------------------------*\
                       Class shiftedForce Declaration
\*---------------------------------------------------------------------------*/

class shiftedForce
:
    public energyScalingFunction
{
    // Private data

        scalar rCut_;

        scalar e_at_rCut_;

        scalar de_dr_at_rCut_;


public:

    //- Runtime type information
    TypeName("shiftedForce");


    // Constructors

        //- Construct from components
        shiftedForce
        (
            const word& name,
            const dictionary& energyScalingFunctionProperties,
            const pairPotential& pairPot
        );


    // Destructor

        ~shiftedForce()
        {}


    // Member Functions

        void scaleEnergy(scalar& e, const scalar r) const;

        //- Read dictionary
        bool read(const dictionary& energyScalingFunctionProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace energyScalingFunctions
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
