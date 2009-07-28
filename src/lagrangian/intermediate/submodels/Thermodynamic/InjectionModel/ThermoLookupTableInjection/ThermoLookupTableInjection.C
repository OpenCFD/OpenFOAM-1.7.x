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

\*---------------------------------------------------------------------------*/

#include "ThermoLookupTableInjection.H"
#include "scalarIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::ThermoLookupTableInjection<CloudType>::INPUT_FILE_COLS = 11;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::ThermoLookupTableInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return round(injectorCells_.size()*(time1 - time0)*nParcelsPerSecond_);
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::ThermoLookupTableInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
) const
{
    scalar volume = 0.0;
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        forAll(mDot_, injectorI)
        {
            volume += mDot_[injectorI]/rho_[injectorI]*(time1 - time0);
        }
    }

    return volume;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoLookupTableInjection<CloudType>::ThermoLookupTableInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    inputFileName_(this->coeffDict().lookup("inputFile")),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    nParcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
    ),
    x_(0),
    U_(0),
    d_(0),
    rho_(0),
    mDot_(0),
    T_(0),
    cp_(0),
    injectorCells_(0)
{
    scalarListIOList injectorData
    (
        IOobject
        (
            inputFileName_,
            owner.db().time().constant(),
            owner.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    x_.setSize(injectorData.size());
    U_.setSize(injectorData.size());
    d_.setSize(injectorData.size());
    rho_.setSize(injectorData.size());
    mDot_.setSize(injectorData.size());
    T_.setSize(injectorData.size());
    cp_.setSize(injectorData.size());

    // Populate lists
    forAll(injectorData, injectorI)
    {
        if (injectorData[injectorI].size() != INPUT_FILE_COLS)
        {
            FatalErrorIn
            (
                "ThermoLookupTableInjection"
                "("
                    "const dictionary&,"
                    "CloudType& owner"
                ")"
            )   << "Incorrect number of entries in injector specification "
                << "- found " << injectorData[injectorI].size()
                << ", expected " << INPUT_FILE_COLS << ":" << nl
                << "    x0 x1 x2 u0 u1 u2 d rho mDot T cp"
                << nl << exit(FatalError);
        }
        x_[injectorI].component(0) = injectorData[injectorI][0];
        x_[injectorI].component(1) = injectorData[injectorI][1];
        x_[injectorI].component(2) = injectorData[injectorI][2];
        U_[injectorI].component(0) = injectorData[injectorI][3];
        U_[injectorI].component(1) = injectorData[injectorI][4];
        U_[injectorI].component(2) = injectorData[injectorI][5];
        d_[injectorI] = injectorData[injectorI][6];
        rho_[injectorI] = injectorData[injectorI][7];
        mDot_[injectorI] = injectorData[injectorI][8];
        T_[injectorI] = injectorData[injectorI][9];
        cp_[injectorI] = injectorData[injectorI][10];
   }

    // Set/cache the injector cells
    injectorCells_.setSize(injectorData.size());
    forAll(x_, injectorI)
    {
        this->findCellAtPosition(injectorCells_[injectorI], x_[injectorI]);
    }

    // Determine volume of particles to inject
    this->volumeTotal_ = 0.0;
    forAll(mDot_, injectorI)
    {
        this->volumeTotal_ += mDot_[injectorI]/rho_[injectorI];
    }
    this->volumeTotal_ *= duration_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoLookupTableInjection<CloudType>::~ThermoLookupTableInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ThermoLookupTableInjection<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::ThermoLookupTableInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
void Foam::ThermoLookupTableInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label nParcels,
    const scalar time,
    vector& position,
    label& cellOwner
)
{
    label injectorI = parcelI*injectorCells_.size()/nParcels;

    position = x_[injectorI];
    cellOwner = injectorCells_[injectorI];
}


template<class CloudType>
void Foam::ThermoLookupTableInjection<CloudType>::setProperties
(
    const label parcelI,
    const label nParcels,
    const scalar,
    typename CloudType::parcelType* pPtr
)
{
    label injectorI = parcelI*injectorCells_.size()/nParcels;

    // set particle velocity
    parcel.U() = U_[injectorI];

    // set particle diameter
    parcel.d() = d_[injectorI];

    // set particle density
    parcel.rho() = rho_[injectorI];

    // set particle temperature
    parcel.T() = T_[injectorI];

    // set particle specific heat capacity
    parcel.cp() = cp_[injectorI];
}


template<class CloudType>
bool Foam::ThermoLookupTableInjection<CloudType>::fullyDescribed() const
{
    return true;
}


template<class CloudType>
bool Foam::ThermoLookupTableInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
