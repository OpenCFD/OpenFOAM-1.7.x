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
    Foam::ReactingCloud

Description
    Templated base class for reacting cloud

    - Adds to thermodynamic cloud
      - Variable composition (single phase)
      - Phase change

SourceFiles
    ReactingCloudI.H
    ReactingCloud.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingCloud_H
#define ReactingCloud_H

#include "ThermoCloud.H"
#include "reactingCloud.H"
#include "multiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class CompositionModel;

template<class CloudType>
class PhaseChangeModel;

/*---------------------------------------------------------------------------*\
                      Class ReactingCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingCloud
:
    public ThermoCloud<ParcelType>,
    public reactingCloud
{
public:

    //- Type of thermodynamics the cloud was instantiated for
    typedef typename ParcelType::thermoType thermoType;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        ReactingCloud(const ReactingCloud&);

        //- Disallow default bitwise assignment
        void operator=(const ReactingCloud&);


protected:

    // Protected data

        //- Parcel constant properties
        typename ParcelType::constantProperties constProps_;

        //- Multi-component carrier phase thermo
        multiComponentMixture<thermoType>& mcCarrierThermo_;


        // References to the cloud sub-models

            //- Reacting composition model
            autoPtr<CompositionModel<ReactingCloud<ParcelType> > >
                compositionModel_;

            //- Reacting phase change model
            autoPtr<PhaseChangeModel<ReactingCloud<ParcelType> > >
                phaseChangeModel_;


        // Sources

            //- Mass transfer fields - one per carrier phase specie
            PtrList<DimensionedField<scalar, volMesh> > rhoTrans_;


        // Check

            //- Total mass transferred to continuous phase via phase change
            scalar dMassPhaseChange_;


    // Protected Member Functions

        // New parcel helper functions

            //- Check that size of a composition field is valid
            void checkSuppliedComposition
            (
                const scalarField& YSupplied,
                const scalarField& Y,
                const word& YName
            );


        // Cloud evolution functions

            //- Pre-evolve
            void preEvolve();

            //- Post-evolve
            void postEvolve();


public:

    // Constructors

        //- Construct given carrier gas fields
        ReactingCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            basicThermo& thermo
        );


    //- Destructor
    virtual ~ReactingCloud();


    //- Type of parcel the cloud was instantiated for
    typedef ParcelType parcelType;


    // Member Functions

        // Access

            //- Return the constant properties
            inline const typename ParcelType::constantProperties&
                constProps() const;

            //- Return const access to multi-component carrier phase thermo
            inline const multiComponentMixture<thermoType>&
                mcCarrierThermo() const;

            //- Return access to multi-component carrier phase thermo
            inline multiComponentMixture<thermoType>& mcCarrierThermo();


            // Sub-models

                //- Return reference to reacting composition model
                inline const CompositionModel<ReactingCloud<ParcelType> >&
                    composition() const;

                //- Return reference to reacting phase change model
                inline const PhaseChangeModel<ReactingCloud<ParcelType> >&
                    phaseChange() const;


            // Sources

                //- Mass

                    //- Return reference to mass source for field i
                    inline DimensionedField<scalar, volMesh>&
                        rhoTrans(const label i);

                    //- Return reference to mass source fields
                    inline PtrList<DimensionedField<scalar, volMesh> >&
                        rhoTrans();

                    //- Return tmp mass source for field i - fully explicit
                    inline tmp<DimensionedField<scalar, volMesh> >
                        Srho(const label i) const;

                    //- Return tmp total mass source for carrier phase
                    //  - fully explicit
                    inline tmp<DimensionedField<scalar, volMesh> > Srho() const;


        // Check

            //- Print cloud information
            void info() const;

            //- Add to cumulative phase change mass transfer
            void addToMassPhaseChange(const scalar dMass);


        // Cloud evolution functions

            //- Check parcel properties
            void checkParcelProperties
            (
                ParcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Reset the spray source terms
            void resetSourceTerms();

            //- Evolve the spray (inject, move)
            void evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReactingCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
