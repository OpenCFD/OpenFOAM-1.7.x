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

\*---------------------------------------------------------------------------*/

#include "HeatTransferModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::HeatTransferModel<CloudType>::HeatTransferModel(CloudType& owner)
:
    dict_(dictionary::null),
    owner_(owner),
    coeffDict_(dictionary::null)
{}


template<class CloudType>
Foam::HeatTransferModel<CloudType>::HeatTransferModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    dict_(dict),
    owner_(owner),
    coeffDict_(dict.subDict(type + "Coeffs"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::HeatTransferModel<CloudType>::~HeatTransferModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
const CloudType& Foam::HeatTransferModel<CloudType>::owner() const
{
    return owner_;
}


template<class CloudType>
const Foam::dictionary& Foam::HeatTransferModel<CloudType>::dict() const
{
    return dict_;
}


template<class CloudType>
const Foam::dictionary& Foam::HeatTransferModel<CloudType>::coeffDict() const
{
    return coeffDict_;
}


template<class CloudType>
Foam::scalar Foam::HeatTransferModel<CloudType>::h
(
    const scalar dp,
    const vector& Ur,
    const scalar rhoc,
    const scalar rhop,
    const scalar cpc,
    const scalar cpp,
    const scalar muc
) const
{
    const scalar Re = rhoc*mag(Ur)*dp/(muc + ROOTVSMALL);

//    const scalar Pr = muc/alphac;
    const scalar Pr = this->Pr();

    const scalar Nu = this->Nu(Re, Pr);

    const scalar kappa = cpc*muc/Pr;

    return Nu*kappa/dp;
}


template<class CloudType>
Foam::scalar Foam::HeatTransferModel<CloudType>::Cu
(
    const scalar dp,
    const vector& Ur,
    const scalar rhoc,
    const scalar rhop,
    const scalar cpc,
    const scalar cpp,
    const scalar muc
) const
{
    const scalar Re = rhoc*mag(Ur)*dp/(muc + ROOTVSMALL);

//    const scalar Pr = muc/alphac;
    const scalar Pr = this->Pr();

    const scalar Nu = this->Nu(Re, Pr);

    const scalar kappa = cpc*muc/Pr;

    const scalar htc = Nu*kappa/dp;

    return 6.0*htc/(dp*rhop*cpp);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewHeatTransferModel.C"


// ************************************************************************* //

