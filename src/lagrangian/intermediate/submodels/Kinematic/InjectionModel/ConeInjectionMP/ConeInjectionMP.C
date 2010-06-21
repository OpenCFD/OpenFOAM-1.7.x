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

\*---------------------------------------------------------------------------*/

#include "ConeInjectionMP.H"
#include "DataEntry.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::ConeInjectionMP<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        const scalar targetVolume = volumeFlowRate_().integrate(0, time1);

        const label targetParcels =
            parcelsPerInjector_*targetVolume/this->volumeTotal_;

        const label nToInject = targetParcels - nInjected_;

        nInjected_ += nToInject;

        return positions_.size()*nToInject;
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::ConeInjectionMP<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return volumeFlowRate_().integrate(time0, time1);
    }
    else
    {
        return 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeInjectionMP<CloudType>::ConeInjectionMP
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    positionsFile_(this->coeffDict().lookup("positionsFile")),
    positions_
    (
        IOobject
        (
            positionsFile_,
            owner.db().time().constant(),
            owner.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    injectorCells_(positions_.size()),
    axesFile_(this->coeffDict().lookup("axesFile")),
    axes_
    (
        IOobject
        (
            axesFile_,
            owner.db().time().constant(),
            owner.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    parcelsPerInjector_
    (
        readScalar(this->coeffDict().lookup("parcelsPerInjector"))
    ),
    volumeFlowRate_
    (
        DataEntry<scalar>::New
        (
            "volumeFlowRate",
            this->coeffDict()
        )
    ),
    Umag_
    (
        DataEntry<scalar>::New
        (
            "Umag",
            this->coeffDict()
        )
    ),
    thetaInner_
    (
        DataEntry<scalar>::New
        (
            "thetaInner",
            this->coeffDict()
        )
    ),
    thetaOuter_
    (
        DataEntry<scalar>::New
        (
            "thetaOuter",
            this->coeffDict()
        )
    ),
    parcelPDF_
    (
        pdfs::pdf::New
        (
            this->coeffDict().subDict("parcelPDF"),
            owner.rndGen()
        )
    ),
    nInjected_(this->parcelsAddedTotal()),
    tanVec1_(positions_.size()),
    tanVec2_(positions_.size())
{
    // Normalise direction vector and determine direction vectors
    // tangential to direction
    forAll(axes_, i)
    {
        axes_[i] /= mag(axes_[i]);

        vector tangent = vector::zero;
        scalar magTangent = 0.0;

        while (magTangent < SMALL)
        {
            vector v = this->owner().rndGen().vector01();

            tangent = v - (v & axes_[i])*axes_[i];
            magTangent = mag(tangent);
        }

        tanVec1_[i] = tangent/magTangent;
        tanVec2_[i] = axes_[i]^tanVec1_[i];
    }

    // Set total volume to inject
    this->volumeTotal_ = volumeFlowRate_().integrate(0.0, duration_);

    // Set/cache the injector cells
    forAll(positions_, i)
    {
        this->findCellAtPosition
        (
            injectorCells_[i],
            positions_[i]
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ConeInjectionMP<CloudType>::~ConeInjectionMP()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ConeInjectionMP<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::ConeInjectionMP<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
void Foam::ConeInjectionMP<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label,
    const scalar,
    vector& position,
    label& cellOwner
)
{
    const label i = parcelI%positions_.size();

    position = positions_[i];
    cellOwner = injectorCells_[i];
}


template<class CloudType>
void Foam::ConeInjectionMP<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar time,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    const label i = parcelI%positions_.size();

    const scalar deg2Rad = mathematicalConstant::pi/180.0;

    scalar t = time - this->SOI_;
    scalar ti = thetaInner_().value(t);
    scalar to = thetaOuter_().value(t);
    scalar coneAngle = this->owner().rndGen().scalar01()*(to - ti) + ti;

    coneAngle *= deg2Rad;
    scalar alpha = sin(coneAngle);
    scalar dcorr = cos(coneAngle);
    scalar beta = mathematicalConstant::twoPi*this->owner().rndGen().scalar01();

    vector normal = alpha*(tanVec1_[i]*cos(beta) + tanVec2_[i]*sin(beta));
    vector dirVec = dcorr*axes_[i];
    dirVec += normal;

    dirVec /= mag(dirVec);

    parcel.U() = Umag_().value(t)*dirVec;

    // set particle diameter
    parcel.d() = parcelPDF_().sample();
}


template<class CloudType>
bool Foam::ConeInjectionMP<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::ConeInjectionMP<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
