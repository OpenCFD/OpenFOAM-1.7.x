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

\*---------------------------------------------------------------------------*/

#include "PatchInjection.H"
#include "DataEntry.H"
#include "pdf.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::PatchInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return round(fraction_*(time1 - time0)*parcelsPerSecond_);
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::PatchInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return fraction_*volumeFlowRate_().integrate(time0, time1);
    }
    else
    {
        return 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchInjection<CloudType>::PatchInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    patchName_(this->coeffDict().lookup("patchName")),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    parcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
    ),
    U0_(this->coeffDict().lookup("U0")),
    volumeFlowRate_
    (
        DataEntry<scalar>::New
        (
            "volumeFlowRate",
            this->coeffDict()
        )
    ),
    parcelPDF_
    (
        pdf::New
        (
            this->coeffDict().subDict("parcelPDF"),
            owner.rndGen()
        )
    ),
    cellOwners_(),
    fraction_(1.0)
{
    label patchId = owner.mesh().boundaryMesh().findPatchID(patchName_);

    if (patchId < 0)
    {
        FatalErrorIn
        (
            "PatchInjection<CloudType>::PatchInjection"
            "("
                "const dictionary&, "
                "CloudType&"
            ")"
        )   << "Requested patch " << patchName_ << " not found" << nl
            << "Available patches are: " << owner.mesh().boundaryMesh().names()
            << nl << exit(FatalError);
    }

    const polyPatch& patch = owner.mesh().boundaryMesh()[patchId];

    cellOwners_ = patch.faceCells();

    label patchSize = cellOwners_.size();
    label totalPatchSize = patchSize;
    reduce(totalPatchSize, sumOp<label>());
    fraction_ = scalar(patchSize)/totalPatchSize;

    // Set total volume/mass to inject
    this->volumeTotal_ = fraction_*volumeFlowRate_().integrate(0.0, duration_);
    this->massTotal_ *= fraction_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchInjection<CloudType>::~PatchInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::PatchInjection<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::PatchInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
void Foam::PatchInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector& position,
    label& cellOwner
)
{
    if (cellOwners_.size() > 0)
    {
        label cellI = this->owner().rndGen().integer(0, cellOwners_.size() - 1);

        cellOwner = cellOwners_[cellI];
        position = this->owner().mesh().C()[cellOwner];
    }
    else
    {
        cellOwner = -1;
        // dummy position
        position = pTraits<vector>::max;
    }
}


template<class CloudType>
void Foam::PatchInjection<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    parcel.U() = U0_;

    // set particle diameter
    parcel.d() = parcelPDF_->sample();
}


template<class CloudType>
bool Foam::PatchInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::PatchInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
