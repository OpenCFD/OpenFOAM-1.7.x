/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    coeffDict_(dictionary::null),
    BirdCorrection_(false)
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
    coeffDict_(dict.subDict(type + "Coeffs")),
    BirdCorrection_(coeffDict_.lookup("BirdCorrection"))
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
const Foam::Switch& Foam::HeatTransferModel<CloudType>::BirdCorrection() const
{
    return BirdCorrection_;
}


template<class CloudType>
Foam::scalar Foam::HeatTransferModel<CloudType>::htc
(
    const scalar dp,
    const scalar Re,
    const scalar Pr,
    const scalar kappa,
    const scalar NCpW
) const
{
    const scalar Nu = this->Nu(Re, Pr);

    scalar htc = Nu*kappa/dp;

    if (BirdCorrection_ && (mag(htc) > ROOTVSMALL) && (mag(NCpW) > ROOTVSMALL))
    {
        const scalar phit = min(NCpW/htc, 50);
        if (phit > 0.001)
        {
            htc *= phit/(exp(phit) - 1.0);
        }
    }

    return htc;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "NewHeatTransferModel.C"


// ************************************************************************* //

