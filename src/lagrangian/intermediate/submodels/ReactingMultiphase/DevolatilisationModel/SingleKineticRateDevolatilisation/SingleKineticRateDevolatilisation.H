/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::SingleKineticRateDevolatilisation

Description
    Single kinetic rate devolatisation model

\*---------------------------------------------------------------------------*/

#ifndef SingleKineticRateDevolatilisation_H
#define SingleKineticRateDevolatilisation_H

#include "DevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
              Class SingleKineticRateDevolatilisation Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SingleKineticRateDevolatilisation
:
    public DevolatilisationModel<CloudType>
{
    // Private data

        // Model constants

            //- Activation energy
            const scalar A1_;

            //- Pre-exoponential factor
            const scalar E_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar volatileResidualCoeff_;


public:

    //- Runtime type information
    TypeName("SingleKineticRateDevolatilisation");


    // Constructors

        //- Construct from dictionary
        SingleKineticRateDevolatilisation
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~SingleKineticRateDevolatilisation();


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
#   include "SingleKineticRateDevolatilisation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
