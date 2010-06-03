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

#include "FieldActivatedInjection.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::FieldActivatedInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if (sum(nParcelsInjected_) < nParcelsPerInjector_*positions_.size())
    {
        return positions_.size();
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::FieldActivatedInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if (sum(nParcelsInjected_) < nParcelsPerInjector_*positions_.size())
    {
        return this->volumeTotal_/nParcelsPerInjector_;
    }
    else
    {
        return 0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FieldActivatedInjection<CloudType>::FieldActivatedInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    factor_(readScalar(this->coeffDict().lookup("factor"))),
    referenceField_
    (
        owner.db().objectRegistry::lookupObject<volScalarField>
        (
            this->coeffDict().lookup("referenceField")
        )
    ),
    thresholdField_
    (
        owner.db().objectRegistry::lookupObject<volScalarField>
        (
            this->coeffDict().lookup("thresholdField")
        )
    ),
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
    nParcelsPerInjector_
    (
        readLabel(this->coeffDict().lookup("parcelsPerInjector"))
    ),
    nParcelsInjected_(positions_.size(), 0),
    U0_(this->coeffDict().lookup("U0")),
    diameters_(positions_.size()),
    parcelPDF_
    (
        pdf::New
        (
            this->coeffDict().subDict("parcelPDF"),
            owner.rndGen()
        )
    )
{
    // Construct parcel diameters - one per injector cell
    forAll(diameters_, i)
    {
        diameters_[i] = parcelPDF_->sample();
    }

    // Determine total volume of particles to inject
    this->volumeTotal_ =
        nParcelsPerInjector_*sum(pow3(diameters_))*mathematicalConstant::pi/6.0;

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
Foam::FieldActivatedInjection<CloudType>::~FieldActivatedInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::FieldActivatedInjection<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::FieldActivatedInjection<CloudType>::timeEnd() const
{
    return GREAT;
}


template<class CloudType>
void Foam::FieldActivatedInjection<CloudType>::setPositionAndCell
(
    const label parcelI,
    const label,
    const scalar,
    vector& position,
    label& cellOwner
)
{
    position = positions_[parcelI];
    cellOwner = injectorCells_[parcelI];
}


template<class CloudType>
void Foam::FieldActivatedInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    parcel.U() = U0_;

    // set particle diameter
    parcel.d() = diameters_[parcelI];
}


template<class CloudType>
bool Foam::FieldActivatedInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::FieldActivatedInjection<CloudType>::validInjection
(
    const label parcelI
)
{
    const label cellI = injectorCells_[parcelI];

    if
    (
        nParcelsInjected_[parcelI] < nParcelsPerInjector_
     && factor_*referenceField_[cellI] > thresholdField_[cellI]
    )
    {
        nParcelsInjected_[parcelI]++;
        return true;
    }

    return false;
}


// ************************************************************************* //
