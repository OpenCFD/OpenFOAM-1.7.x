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
    Foam::spray

Description
    A spray is a cloud of parcels

\*---------------------------------------------------------------------------*/

#ifndef spray_H
#define spray_H

#include "parcel.H"
#include "injector.H"
#include "IOPtrList.H"
#include "interpolation.H"
#include "liquid.H"
#include "autoPtr.H"
#include "liquidMixture.H"
#include "Random.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class atomizationModel;
class breakupModel;
class collisionModel;
class dispersionModel;
class dragModel;
class evaporationModel;
class injectorModel;
class heatTransferModel;
class wallModel;

class basicMultiComponentMixture;

/*---------------------------------------------------------------------------*\
                        Class spray Declaration
\*---------------------------------------------------------------------------*/

class spray
:
    public Cloud<parcel>
{
    // Private data

        // References to the database and meshes

            const Time& runTime_;
            scalar time0_;
            const fvMesh& mesh_;

            //- Random number generator
            Random rndGen_;

        //- Acceleration due to gravity
        const vector& g_;


        // References to the physical fields

            const volVectorField& U_;
            const volScalarField& rho_;
            const volScalarField& p_;
            const volScalarField& T_;


        //- The spray properties
        IOdictionary sprayProperties_;


        //- Ambient Pressure
        scalar ambientPressure_;

        //- Ambient Temperature
        scalar ambientTemperature_;


        //- The injectors
        IOPtrList<injector> injectors_;


        // References to the spray sub-models

            autoPtr<atomizationModel> atomization_;
            autoPtr<dragModel> drag_;
            autoPtr<evaporationModel> evaporation_;
            autoPtr<heatTransferModel> heatTransfer_;
            autoPtr<wallModel> wall_;
            autoPtr<breakupModel> breakupModel_;
            autoPtr<collisionModel> collisionModel_;
            autoPtr<dispersionModel> dispersionModel_;
            autoPtr<liquidMixture> fuels_;
            autoPtr<injectorModel> injectorModel_;


        //- Minimum number of lagrangian subcycles
        const label subCycles_;


        // Composition properties

            const PtrList<gasThermoPhysics>& gasProperties_;
            const basicMultiComponentMixture& composition_;

            List<label> liquidToGasIndex_;
            List<label> gasToLiquidIndex_;
            List<bool> isLiquidFuel_;


        // Necessary 2D-information

            bool twoD_;
            vector axisOfSymmetry_;
            vector axisOfWedge_;
            vector axisOfWedgeNormal_;
            scalar angleOfWedge_;


        // Interpolation

            dictionary interpolationSchemes_;

            autoPtr<interpolation<vector> > UInterpolator_;
            autoPtr<interpolation<scalar> > rhoInterpolator_;
            autoPtr<interpolation<scalar> > pInterpolator_;
            autoPtr<interpolation<scalar> > TInterpolator_;


        // Spray Source Terms

            //- Momentum
            vectorField sms_;

            //- Enthalpy
            scalarField shs_;

            //- Mass
            PtrList<scalarField> srhos_;

            //- The total mass of the injected liquid
            scalar totalInjectedLiquidMass_;

            //- The (total added) injected kinetic energy of the liquid
            scalar injectedLiquidKE_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        spray(const spray&);

        //- Disallow default bitwise assignment
        void operator=(const spray&);


public:

    // Constructors

        //- Construct from components
        spray
        (
            const volVectorField& U,
            const volScalarField& rho,
            const volScalarField& p,
            const volScalarField& T,
            const basicMultiComponentMixture& composition,
            const PtrList<gasThermoPhysics>& gasProperties,
            const dictionary& thermophysicalProperties,
            const dimensionedVector& g,
            bool readFields = true
        );


    // Destructor

        ~spray();


    // Member Functions

        // Spray tracking and evolution functions

            //- Evolve the spray (move, inject and breakup)
            void evolve();

            //- Move the spray parcels
            void move();

            //- Inject more parcels
            void inject();

            //- Primary breakup droplets
            void atomizationLoop();


            //- Secondary breakup droplets
            void breakupLoop();


        // Access

            inline const Time& runTime() const;
            inline const fvMesh& mesh() const;

            inline const volVectorField& U() const;
            inline const volScalarField& rho() const;
            inline const volScalarField& p() const;
            inline const volScalarField& T() const;

            inline PtrList<injector>& injectors();
            inline const PtrList<injector>& injectors() const;

            inline const atomizationModel& atomization() const;
            inline const breakupModel& breakup() const;
            inline const collisionModel& collisions() const;
            inline const dispersionModel& dispersion() const;
            inline const dragModel& drag() const;
            inline const evaporationModel& evaporation() const;
            inline const heatTransferModel& heatTransfer() const;
            inline const injectorModel& injection() const;
            inline const wallModel& wall() const;

            inline tmp<volVectorField> momentumSource() const;
            inline tmp<volScalarField> evaporationSource(const label i) const;
            inline tmp<volScalarField> heatTransferSource() const;

            inline Random& rndGen();
            inline label subCycles() const;
            inline const vector& g() const;

            inline const liquidMixture& fuels() const;
            inline const PtrList<gasThermoPhysics>& gasProperties() const;
            inline const basicMultiComponentMixture& composition() const;

            inline const List<label>& liquidToGasIndex() const;
            inline const List<label>& gasToLiquidIndex() const;
            inline const List<bool>& isLiquidFuel() const;

            inline const bool& twoD() const;
            inline const vector& axisOfSymmetry() const;
            inline const vector& axisOfWedge() const;
            inline const vector& axisOfWedgeNormal() const;
            inline const scalar& angleOfWedge() const;

            inline const interpolation<vector>& UInterpolator() const;
            inline const interpolation<scalar>& rhoInterpolator() const;
            inline const interpolation<scalar>& pInterpolator() const;
            inline const interpolation<scalar>& TInterpolator() const;

            inline vectorField& sms();
            inline const vectorField& sms() const;

            inline scalarField& shs();
            inline const scalarField& shs() const;

            inline PtrList<scalarField>& srhos();
            inline const PtrList<scalarField>& srhos() const;

            inline const scalar& ambientPressure() const;

            inline const scalar& ambientTemperature() const;


        // Check

            //- Returns the liquid mass that has been injected
            scalar injectedMass(const scalar t) const;

           //- Returns the liquid mass that will be injected by the injectors
            scalar totalMassToInject() const;

            //- Returns the injected enthalpy
            scalar injectedEnthalpy(const scalar t) const;

            //- Returns current total liquid mass in the domain
            scalar liquidMass() const;

            //- Returns the enthalpy of all the liquid in the domain
            // Hdrop = Hgas - Hlat
            scalar liquidEnthalpy() const;

            //- Returns the enthalpy (total) of all the liquid in the domain
            // Hdrop = Hgas - Hlat + (P-Psat)/rhoDrop;
            scalar liquidTotalEnthalpy() const;

            //- Returns the kinetic energy of the liquid phase
            scalar liquidKineticEnergy() const;

            //- Returns the injected kinetic energy of the liquid phase
            scalar injectedLiquidKineticEnergy() const;

            //- Returns the droplet penetration for 'prc' percent of the
            //  liquid from nozzle 'nozzlei'
            scalar liquidPenetration
            (
                const label nozzlei,
                const scalar prc
            ) const;

            //- Returns the droplet penetration for 'prc' percent of the
            //  liquid from nozzle 0
            scalar liquidPenetration(const scalar prc) const;

            //- Return Sauter Mean Diameter
            scalar smd() const;

            //- Return Maximum Diameter
            scalar maxD() const;

            //- Return Ambient Pressure
            void calculateAmbientPressure();

            //- Return Ambient Temperature
            void calculateAmbientTemperature();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sprayI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
