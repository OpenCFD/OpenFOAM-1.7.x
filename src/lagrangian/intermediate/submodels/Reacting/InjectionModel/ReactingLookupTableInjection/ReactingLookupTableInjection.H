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
    Foam::ReactingLookupTableInjection

Description
    Particle injection sources read from look-up table. Each row corresponds to
    an injection site.

    (
        (x y z u v w d rho mDot T cp Y0..YN)
    );

    where:
        x, y, z = global cartesian co-ordinates [m]
        u, v, w = global cartesian velocity components [m/s]
        d       = diameter [m]
        rho     = density [kg/m3]
        mDot    = mass flow rate [kg/m3]
        T       = temperature [K]
        cp      = specific heat capacity [J/kg/K]

SourceFiles
    ReactingLookupTableInjection.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingLookupTableInjection_H
#define ReactingLookupTableInjection_H

#include "InjectionModel.H"
#include "pdf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class ReactingLookupTableInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ReactingLookupTableInjection
:
    public InjectionModel<CloudType>
{
    // Private data

        //- Name of file containing injector/parcel data
        const word inputFileName_;

        //- Injection duration - common to all injection sources
        const scalar duration_;

        //- Number of parcels per injector - common to all injection sources
        const label nParcelsPerSecond_;

        //- List of parcel position per injector / [m]
        List<point> x_;

        //- List of parcel velocity per injector / [m]
        List<vector> U_;

        //- List of parcel diameter per injector / [m]
        List<scalar> d_;

        //- List of parcel fluid density pre injector / [kg/m3]
        List<scalar> rho_;

        //- List of parcel injection mass flow per injector / [kg/s]
        List<scalar> mDot_;

        //- List of parcel temperature flow per injector / [K]
        List<scalar> T_;

        //- List of parcel specific heat capacity per injector / [J/kg/K]
        List<scalar> cp_;

        //- List of parcel specie mass fractions per injector / []
        List<scalarList> Y_;

        //- List of injector cells per injector
        List<label> injectorCells_;

        //- Number of columns expected in input file
        static label INPUT_FILE_COLS;


protected:

    // Protected member functions

        //- Number of parcels to introduce over the time step relative to SOI
        label parcelsToInject
        (
            const scalar time0,
            const scalar time1
        ) const;

        //- Volume of parcels to introduce over the time step relative to SOI
        scalar volumeToInject
        (
            const scalar time0,
            const scalar time1
        ) const;


public:

    //- Runtime type information
    TypeName("ReactingLookupTableInjection");


    // Constructors

        //- Construct from dictionary
        ReactingLookupTableInjection
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~ReactingLookupTableInjection();


    // Member Functions

        //- Flag to indicate whether model activates injection model
        bool active() const;

        //- Return the end-of-injection time
        scalar timeEnd() const;


        // Injection geometry

            //- Set the injection position and owner cell
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReactingLookupTableInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
