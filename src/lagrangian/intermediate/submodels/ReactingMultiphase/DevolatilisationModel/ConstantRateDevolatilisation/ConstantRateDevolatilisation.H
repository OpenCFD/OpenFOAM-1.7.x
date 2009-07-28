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
    Foam::ConstantRateDevolatilisation

Description
    Constant rate devolatisation model
    - need to set vapourisation temperature to 600 K

\*---------------------------------------------------------------------------*/

#ifndef ConstantRateDevolatilisation_H
#define ConstantRateDevolatilisation_H

#include "DevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                 Class ConstantRateDevolatilisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ConstantRateDevolatilisation
:
    public DevolatilisationModel<CloudType>
{
    // Private data

        // Model constants

            //- Rate constant (suggested default = 12) [1/s]
            const scalar A0_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar volatileResidualCoeff_;


public:

    //- Runtime type information
    TypeName("ConstantRateDevolatilisation");


    // Constructors

        //- Construct from dictionary
        ConstantRateDevolatilisation
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~ConstantRateDevolatilisation();


    // Member Functions

        //- Flag to indicate whether model activates devolatilisation model
        virtual bool active() const;

        //- Update model
        virtual scalar calculate
        (
            const scalar dt,
            const scalar mass0,
            const scalar mass,
            const scalar T,
            const scalar YVolatile0,
            const scalar YVolatile,
            bool& canCombust
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ConstantRateDevolatilisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
