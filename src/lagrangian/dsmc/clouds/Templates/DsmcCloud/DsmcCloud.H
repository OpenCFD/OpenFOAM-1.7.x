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
    Foam::DsmcCloud

Description
    Templated base class for dsmc cloud

SourceFiles
    DsmcCloudI.H
    DsmcCloud.C

\*---------------------------------------------------------------------------*/

#ifndef DsmcCloud_H
#define DsmcCloud_H

#include "Cloud.H"
#include "DsmcBaseCloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
#include "fvMesh.H"
#include "volFields.H"
#include "scalarIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class BinaryCollisionModel;

template<class CloudType>
class WallInteractionModel;

template<class CloudType>
class InflowBoundaryModel;

/*---------------------------------------------------------------------------*\
                         Class DsmcCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class DsmcCloud
:
    public Cloud<ParcelType>,
    public DsmcBaseCloud
{
    // Private data

        //- Cloud type - used to set the name of the parcel properties
        //  dictionary by appending "Properties"
        const word cloudName_;

        //- References to the mesh and time databases
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- A list of unique instances of molecule types in the simulation.
        // The position of an entry in the list maps to the label identifying
        // the typeId, i.e. where typeIdList_ = (N2 O2 CO2)
        // N2 has typeId label = 0, O2 = 1, CO2 = 2.
        List<word> typeIdList_;

        //- Number of real atoms/molecules represented by a parcel
        scalar nParticle_;

        //- A data structure holding which particles are in which cell
        List<DynamicList<ParcelType*> > cellOccupancy_;

        //- An IOField holding the value of (sigmaT * cR)max for each cell (see
        // Bird p220). Initialised with the parcels, updated as required, and
        // read in on start/restart.
        volScalarField sigmaTcRMax_;

        //- A field holding the remainder from the previous collision selections
        scalarField collisionSelectionRemainder_;

        //- Heat flux at surface field
        volScalarField q_;

        //- Force density at surface field
        volVectorField fD_;

        //- Parcel constant properties - one for each type
        List<typename ParcelType::constantProperties> constProps_;

        //- Random number generator
        Random rndGen_;

        //- In-cloud cache of deltaT, lookup in submodels and parcel is
        //  expensive
        scalar cachedDeltaT_;


        // References to the macroscopic fields

            //- Temperature
            const volScalarField& T_;

            //- Velocity
            const volVectorField& U_;


        // References to the cloud sub-models

            //- Binary collision model
            autoPtr<BinaryCollisionModel<DsmcCloud<ParcelType> > >
                binaryCollisionModel_;

            //- Wall interaction model
            autoPtr<WallInteractionModel<DsmcCloud<ParcelType> > >
                wallInteractionModel_;

            //- Inflow boundary model
            autoPtr<InflowBoundaryModel<DsmcCloud<ParcelType> > >
                inflowBoundaryModel_;


    // Private Member Functions

        //- Build the constant properties for all of the species
        void buildConstProps();

        //- Record which particles are in which cell
        void buildCellOccupancy();

        //- Initialise the system
        void initialise(const IOdictionary& dsmcInitialiseDict);

        //- Calculate collisions between molecules
        void collisions();

        //- Reset the surface data accumulation field values
        void resetSurfaceDataFields();

        //- Disallow default bitwise copy construct
        DsmcCloud(const DsmcCloud&);

        //- Disallow default bitwise assignment
        void operator=(const DsmcCloud&);


public:

    // Static data members

        //- Boltzmann constant
        static scalar kb;


    // Constructors

        //- Construct given name and mesh, will read Parcels from file
        DsmcCloud
        (
            const word& cloudName,
            const volScalarField& T,
            const volVectorField& U
        );

        //- Construct given name and mesh.  Used to initialise.
        DsmcCloud
        (
            const word& cloudName,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~DsmcCloud();


    // Member Functions

        // Access

            // References to the mesh and databases

                //- Return the cloud type
                inline const word& cloudName() const;

                //- Return refernce to the mesh
                inline const fvMesh& mesh() const;


            // References to the dsmc specific data

                //- Return particle properties dictionary
                inline const IOdictionary& particleProperties() const;

                //- Return the idList
                inline const List<word>& typeIdList() const;

                //- Return the number of real particles represented by one
                //  parcel
                inline scalar nParticle() const;

                //- Return the cell occupancy addressing
                inline const List<DynamicList<ParcelType*> >&
                    cellOccupancy() const;

                //- Return the sigmaTcRMax field.  non-const access to allow
                // updating.
                inline volScalarField& sigmaTcRMax();

                //- Return the collision selection remainder field.  non-const
                // access to allow updating.
                inline scalarField& collisionSelectionRemainder();

                //- Return all of the constant properties
                inline const List<typename ParcelType::constantProperties>&
                    constProps() const;

                //- Return the constant properties of the given typeId
                inline const typename ParcelType::constantProperties&
                    constProps(label typeId) const;

                //- Return refernce to the random object
                inline Random& rndGen();

                //- Store (cache) the current value of deltaT
                inline void storeDeltaT();

                //- Return the cached value of deltaT
                inline scalar cachedDeltaT() const;


            // References to the surface data collection fields

                //- Return heat flux at surface field
                inline const volScalarField& q() const;

                //- Return non-const heat flux at surface field
                inline volScalarField& q();

                //- Return force density at surface field
                inline const volVectorField& fD() const;

                //- Return non-const force density at surface field
                inline volVectorField& fD();


            // References to the macroscopic fields

                //- Return macroscopic temperature
                inline const volScalarField& T() const;

                //- Return macroscopic velocity
                inline const volVectorField& U() const;


            // Kinetic theory helper functions

                //- Generate a random velocity sampled from the Maxwellian speed
                // distribution
                vector equipartitionLinearVelocity
                (
                    scalar temperature,
                    scalar mass
                );

                //- Generate a random internal energy, sampled from the
                // equilibrium distribution (Bird eqn 11.22 and 11.23 and
                // adapting code from DSMC3.FOR)
                scalar equipartitionInternalEnergy
                (
                    scalar temperature,
                    scalar internalDegreesOfFreedom
                );


                // From the Maxwellian distribution:
                //- Average particle speed
                inline scalar maxwellianAverageSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                inline scalarField maxwellianAverageSpeed
                (
                    scalarField temperature,
                    scalar mass
                ) const;

                //- RMS particle speed
                inline scalar maxwellianRMSSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                inline scalarField maxwellianRMSSpeed
                (
                    scalarField temperature,
                    scalar mass
                ) const;

                //- Most probable speed
                inline scalar maxwellianMostProbableSpeed
                (
                    scalar temperature,
                    scalar mass
                ) const;

                inline scalarField maxwellianMostProbableSpeed
                (
                    scalarField temperature,
                    scalar mass
                ) const;

            // Sub-models

                //- Return reference to binary elastic collision model
                inline const BinaryCollisionModel<DsmcCloud<ParcelType> >&
                    binaryCollision() const;

                //- Return non-const reference to binary elastic collision model
                inline BinaryCollisionModel<DsmcCloud<ParcelType> >&
                    binaryCollision();

                //- Return reference to wall interaction model
                inline const WallInteractionModel<DsmcCloud<ParcelType> >&
                    wallInteraction() const;

                //- Return non-const reference to wall interaction model
                inline WallInteractionModel<DsmcCloud<ParcelType> >&
                    wallInteraction();

                //- Return reference to wall interaction model
                inline const InflowBoundaryModel<DsmcCloud<ParcelType> >&
                    inflowBoundary() const;

                //- Return non-const reference to wall interaction model
                inline InflowBoundaryModel<DsmcCloud<ParcelType> >&
                    inflowBoundary();


        // Check

            //- Total mass injected
            inline scalar massInjected() const;

            //- Total mass in system
            inline scalar massInSystem() const;

            //- Total linear momentum of the system
            inline vector linearMomentumOfSystem() const;

            //- Total linear kinetic energy in the system
            inline scalar linearKineticEnergyOfSystem() const;

            //- Total internal energy in the system
            inline scalar internalEnergyOfSystem() const;

            //- Print cloud information
            void info() const;

            //- Dump particle positions to .obj file
            void dumpParticlePositions() const;


            // Fields

                //- Return the real particle number density field
                inline const tmp<volScalarField> rhoN() const;

                //- Return the particle mass density field
                inline const tmp<volScalarField> rhoM() const;

                //- Return the field of number of DSMC particles
                inline const tmp<volScalarField> dsmcRhoN() const;

                //- Return the momentum density field
                inline const tmp<volVectorField> momentum() const;

                //- Return the total linear kinetic energy (translational and
                // thermal density field
                inline const tmp<volScalarField> linearKE() const;

                //- Return the internal energy density field
                inline const tmp<volScalarField> internalE() const;

                //- Return the average internal degrees of freedom  field
                inline const tmp<volScalarField> iDof() const;


        // Cloud evolution functions

            //- Add new parcel
            void addNewParcel
            (
                const vector& position,
                const vector& U,
                const scalar Ei,
                const label cellId,
                const label typeId
            );

            //- Evolve the cloud (move, collide)
            void evolve();


            //- Clear the Cloud
            inline void clear();


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DsmcCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DsmcCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
