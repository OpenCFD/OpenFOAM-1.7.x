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
    Foam::ReactingMultiphaseLookupTableInjection

Description
    Particle injection sources read from look-up table. Each row corresponds to
    an injection site.

    (
       (x y z) (u v w) d rho mDot T cp (Y0..Y2) (Yg0..YgN) (Yl0..YlN) (Ys0..YsN)
       (x y z) (u v w) d rho mDot T cp (Y0..Y2) (Yg0..YgN) (Yl0..YlN) (Ys0..YsN)
       ...
       (x y z) (u v w) d rho mDot T cp (Y0..Y2) (Yg0..YgN) (Yl0..YlN) (Ys0..YsN)
    );

    where:
        x, y, z  = global cartesian co-ordinates [m]
        u, v, w  = global cartesian velocity components [m/s]
        d        = diameter [m]
        rho      = density [kg/m3]
        mDot     = mass flow rate [kg/m3]
        T        = temperature [K]
        cp       = specific heat capacity [J/kg/K]
        Y(3)     = total mass fraction of gas (Y0), liquid (Y1), solid (Y3)
        Yg(Ngas) = mass fractions of gaseous components
        Yl(Nliq) = mass fractions of liquid components
        Ys(Nsld) = mass fractions of solid components

SourceFiles
    ReactingMultiphaseLookupTableInjection.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingMultiphaseLookupTableInjection_H
#define ReactingMultiphaseLookupTableInjection_H

#include "InjectionModel.H"
#include "reactingMultiphaseParcelInjectionDataIOList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class ReactingMultiphaseLookupTableInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ReactingMultiphaseLookupTableInjection
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

        //- List of injectors
        reactingMultiphaseParcelInjectionDataIOList injectors_;

        //- List of injector cells per injector
        List<label> injectorCells_;


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
    TypeName("ReactingMultiphaseLookupTableInjection");


    // Constructors

        //- Construct from dictionary
        ReactingMultiphaseLookupTableInjection
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~ReactingMultiphaseLookupTableInjection();


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
#   include "ReactingMultiphaseLookupTableInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
