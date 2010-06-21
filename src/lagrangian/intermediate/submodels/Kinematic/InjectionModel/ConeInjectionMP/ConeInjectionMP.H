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
    Foam::ConeInjectionMP

Description
    Cone injection multi-point

    - User specifies
      - time of start of injection
      - injector positions
      - directions (along injection axes)
      - parcel flow rate
      - parcel velocities
      - inner and outer cone angles
    - Parcel diameters obtained by PDF model

SourceFiles
    ConeInjectionMP.C

\*---------------------------------------------------------------------------*/

#ifndef ConeInjectionMP_H
#define ConeInjectionMP_H

#include "InjectionModel.H"
#include "pdf.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class Type>
class DataEntry;

/*---------------------------------------------------------------------------*\
                       Class ConeInjectionMP Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ConeInjectionMP
:
    public InjectionModel<CloudType>
{
    // Private data

        //- Name of file containing positions data
        const word positionsFile_;

        //- Field of injector positions
        vectorIOField positions_;

        //- Field of cell labels corresoponding to injector positions
        labelList injectorCells_;

        //- Name of file containing axes data
        const word axesFile_;

        //- Field of injector positions
        vectorIOField axes_;

        //- Injection duration [s]
        const scalar duration_;

        //- Number of parcels to introduce per injector
        const label parcelsPerInjector_;

        //- Volume flow rate of parcels to introduce relative to SOI [m^3]
        const autoPtr<DataEntry<scalar> > volumeFlowRate_;

        //- Parcel velocity magnitude relative to SOI [m/s]
        const autoPtr<DataEntry<scalar> > Umag_;

        //- Inner cone angle relative to SOI [deg]
        const autoPtr<DataEntry<scalar> > thetaInner_;

        //- Outer cone angle relative to SOI [deg]
        const autoPtr<DataEntry<scalar> > thetaOuter_;

        //- Parcel size PDF model
        const autoPtr<pdfs::pdf> parcelPDF_;

        //- Number of parcels per injector already injected
        mutable label nInjected_;


        // Tangential vectors to the direction vector

            //- First tangential vector
            vectorList tanVec1_;

            //- Second tangential vector
            vectorList tanVec2_;


protected:

    // Protected member functions

        //- Number of parcels to introduce over the time step
        label parcelsToInject
        (
            const scalar time0,
            const scalar time1
        ) const;

        //- Number of parcels to introduce over the time step
        scalar volumeToInject
        (
            const scalar time0,
            const scalar time1
        ) const;


public:

    //- Runtime type information
    TypeName("ConeInjectionMP");


    // Constructors

        //- Construct from dictionary
        ConeInjectionMP
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~ConeInjectionMP();


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
#   include "ConeInjectionMP.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
