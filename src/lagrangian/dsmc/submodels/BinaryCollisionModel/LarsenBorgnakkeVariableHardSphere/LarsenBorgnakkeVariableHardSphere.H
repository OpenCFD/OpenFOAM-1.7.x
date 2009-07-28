/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    Foam::LarsenBorgnakkeVariableHardSphere

Description
    Variable Hard Sphere BinaryCollision Model with Larsen Borgnakke internal
    energy redistribution.  Based on the INELRS subroutine in Bird's DSMC0R.FOR

\*---------------------------------------------------------------------------*/

#ifndef LarsenBorgnakkeVariableHardSphere_H
#define LarsenBorgnakkeVariableHardSphere_H

#include "BinaryCollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class LarsenBorgnakkeVariableHardSphere Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LarsenBorgnakkeVariableHardSphere
:
    public BinaryCollisionModel<CloudType>
{
    // Private data

        //- Reference temperature
        const scalar Tref_;

        //- Relaxation collision number
        const scalar relaxationCollisionNumber_;


    // Private Member Functions

        //- Calculate the energy ratio for distribution to internal degrees of
        // freedom
        scalar energyRatio
        (
            scalar ChiA,
            scalar ChiB
        );


public:

    //- Runtime type information
    TypeName("LarsenBorgnakkeVariableHardSphere");


    // Constructors

        //- Construct from dictionary
        LarsenBorgnakkeVariableHardSphere
        (
            const dictionary& dict,
            CloudType& cloud
        );


    // Destructor
    virtual ~LarsenBorgnakkeVariableHardSphere();


    // Member Functions

        //- Return the collision cross section * relative velocity product
        virtual scalar sigmaTcR
        (
            label typeIdP,
            label typeIdQ,
            const vector& UP,
            const vector& UQ
        ) const;

        //- Apply collision
        virtual void collide
        (
            label typeIdP,
            label typeIdQ,
            vector& UP,
            vector& UQ,
            scalar& EiP,
            scalar& EiQ
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LarsenBorgnakkeVariableHardSphere.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
