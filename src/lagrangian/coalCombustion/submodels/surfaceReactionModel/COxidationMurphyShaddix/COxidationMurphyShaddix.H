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
    COxidationMurphyShaddix

Description
    Limited to C(s) + O2 -> CO2

\*---------------------------------------------------------------------------*/

#ifndef COxidationMurphyShaddix_H
#define COxidationMurphyShaddix_H

#include "SurfaceReactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
template<class CloudType>
class COxidationMurphyShaddix;

/*---------------------------------------------------------------------------*\
                  Class COxidationMurphyShaddix Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class COxidationMurphyShaddix
:
    public SurfaceReactionModel<CloudType>
{
    // Private data

        // Model constants

            //- Reference diffusion constant
            const scalar D0_;

            //- Reference density for reference diffusion constant
            const scalar rho0_;

            //- Reference temperature for reference diffusion constant
            const scalar T0_;

            //- Exponent for diffusion equation
            const scalar Dn_;

            //- Kinetic rate coefficient
            const scalar A_;

            //- Kinetic rate activation energy
            const scalar E_;

            //- Reaction order
            const scalar n_;

            //- Effective molecular weight of gaseous volatiles
            const scalar WVol_;


        // Static constants

            //- Maximum number of iterations
            static label maxIters_;

            //- Tolerance used in inner iterations
            static scalar tolerance_;


        // Addressing

            //- Cs positions in global/local lists
            label CsLocalId_;

            //- O2 position in global list
            label O2GlobalId_;

            //- CO2 positions in global list
            label CO2GlobalId_;


        // Local copies of thermo properties

            //- Molecular weight of C [kg/kmol]
            scalar WC_;

            //- Molecular weight of O2 [kg/kmol]
            scalar WO2_;


public:

    //- Runtime type information
    TypeName("COxidationMurphyShaddix");


    // Constructors

        //- Construct from dictionary
        COxidationMurphyShaddix
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~COxidationMurphyShaddix();


    // Member Functions

        //- Flag to indicate whether model activates surface reaction model
        virtual bool active() const;

        //- Update surface reactions
        virtual scalar calculate
        (
            const scalar dt,
            const label cellI,
            const scalar d,
            const scalar T,
            const scalar Tc,
            const scalar pc,
            const scalar rhoc,
            const scalar mass,
            const scalarField& YGas,
            const scalarField& YLiquid,
            const scalarField& YSolid,
            const scalarField& YMixture,
            const scalarField& dMassVolatile,
            scalarField& dMassGas,
            scalarField& dMassLiquid,
            scalarField& dMassSolid,
            scalarField& dMassSRCarrier
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "COxidationMurphyShaddix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
