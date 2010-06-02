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
    Foam::kineticTheoryModel

Description

SourceFiles
    kineticTheoryModel.C

\*---------------------------------------------------------------------------*/

#ifndef kineticTheoryModel_H
#define kineticTheoryModel_H

#include "dragModel.H"
#include "phaseModel.H"
#include "autoPtr.H"
#include "viscosityModel.H"
#include "conductivityModel.H"
#include "radialModel.H"
#include "granularPressureModel.H"
#include "frictionalStressModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class kineticTheoryModel Declaration
\*---------------------------------------------------------------------------*/

class kineticTheoryModel
{
    // Private data

        const phaseModel& phasea_;
        const volVectorField& Ua_;
        const volVectorField& Ub_;
        const volScalarField& alpha_;
        const surfaceScalarField& phia_;

        const dragModel& draga_;

        const dimensionedScalar& rhoa_;
        const dimensionedScalar& da_;
        const dimensionedScalar& nua_;

        //- dictionary holding the modeling info
        IOdictionary kineticTheoryProperties_;

        //- use kinetic theory or not.
        Switch kineticTheory_;

        //- use generation == dissipation
        Switch equilibrium_;

        autoPtr<kineticTheoryModels::viscosityModel> viscosityModel_;

        autoPtr<conductivityModel> conductivityModel_;

        autoPtr<radialModel> radialModel_;

        autoPtr<granularPressureModel> granularPressureModel_;

        autoPtr<frictionalStressModel> frictionalStressModel_;

        //- coefficient of restitution
        const dimensionedScalar e_;

        //- maximum packing
        const dimensionedScalar alphaMax_;

        //- min value for which the frictional stresses are zero
        const dimensionedScalar alphaMinFriction_;

        //- material constant for frictional normal stress
        const dimensionedScalar Fr_;

        //- material constant for frictional normal stress
        const dimensionedScalar eta_;

        //- material constant for frictional normal stress
        const dimensionedScalar p_;

        //- angle of internal friction
        const dimensionedScalar phi_;

        //- The granular energy/temperature
        volScalarField Theta_;

        //- The granular viscosity
        volScalarField mua_;

        //- The granular bulk viscosity
        volScalarField lambda_;

        //- The granular pressure
        volScalarField pa_;

        //- The granular temperature conductivity
        volScalarField kappa_;

        //- The radial distribution function
        volScalarField gs0_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        kineticTheoryModel(const kineticTheoryModel&);

        //- Disallow default bitwise assignment
        void operator=(const kineticTheoryModel&);


public:

    // Constructors

        //- Construct from components
        kineticTheoryModel
        (
            const phaseModel& phasea,
            const volVectorField& Ub,
            const volScalarField& alpha,
            const dragModel& draga
        );


    //- Destructor
    virtual ~kineticTheoryModel();


    // Member Functions

        void solve();

        bool on() const
        {
            return kineticTheory_;
        }

        const volScalarField& mua() const
        {
            return mua_;
        }

        const volScalarField& pa() const
        {
            return pa_;
        }

        const volScalarField& lambda() const
        {
            return lambda_;
        }

        const volScalarField& kappa() const
        {
            return kappa_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
