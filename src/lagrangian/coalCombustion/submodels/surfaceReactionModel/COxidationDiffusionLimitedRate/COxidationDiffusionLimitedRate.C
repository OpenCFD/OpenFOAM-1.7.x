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

#include "COxidationDiffusionLimitedRate.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationDiffusionLimitedRate<CloudType>::COxidationDiffusionLimitedRate
(
    const dictionary& dict,
    CloudType& owner
)
:
    SurfaceReactionModel<CloudType>(dict, owner, typeName),
    Sb_(dimensionedScalar(this->coeffDict().lookup("Sb")).value()),
    D_(dimensionedScalar(this->coeffDict().lookup("D")).value()),
    CsLocalId_(-1),
    O2GlobalId_(owner.composition().globalCarrierId("O2")),
    CO2GlobalId_(owner.composition().globalCarrierId("CO2")),
    WC_(0.0),
    WO2_(0.0),
    HcCO2_(0.0)
{
    // Determine Cs ids
    label idSolid = owner.composition().idSolid();
    CsLocalId_ = owner.composition().localId(idSolid, "C");

    // Set local copies of thermo properties
    WO2_ = owner.mcCarrierThermo().speciesData()[O2GlobalId_].W();
    scalar WCO2 = owner.mcCarrierThermo().speciesData()[CO2GlobalId_].W();
    WC_ = WCO2 - WO2_;
    HcCO2_ = owner.mcCarrierThermo().speciesData()[CO2GlobalId_].Hc();

    if (Sb_ < 0)
    {
        FatalErrorIn
        (
            "COxidationDiffusionLimitedRate<CloudType>"
            "("
                "const dictionary&, "
                "CloudType&"
            ")"
        )   << "Stoichiometry of reaction, Sb, must be greater than zero" << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationDiffusionLimitedRate<CloudType>::
~COxidationDiffusionLimitedRate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::COxidationDiffusionLimitedRate<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::COxidationDiffusionLimitedRate<CloudType>::calculate
(
    const scalar dt,
    const label cellI,
    const scalar d,
    const scalar T,
    const scalar Tc,
    const scalar pc,
    const scalar rhoc,
    const scalar mass,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    const scalarField& YMixture,
    const scalarField& dMassVolatile,
    scalarField& dMassGas,
    scalarField& dMassLiquid,
    scalarField& dMassSolid,
    scalarField& dMassSRCarrier
) const
{
    // Fraction of remaining combustible material
    const label idSolid = CloudType::parcelType::SLD;
    const scalar fComb = YMixture[idSolid]*YSolid[CsLocalId_];

    // Surface combustion active combustible fraction is consumed
    if (fComb < SMALL)
    {
        return 0.0;
    }

    // Local mass fraction of O2 in the carrier phase
    const scalar YO2 = this->owner().mcCarrierThermo().Y(O2GlobalId_)[cellI];

    // Change in C mass [kg]
    scalar dmC =
        4.0*mathematicalConstant::pi*d*D_*YO2*Tc*rhoc/(Sb_*(T + Tc))*dt;

    // Limit mass transfer by availability of C
    dmC = min(mass*fComb, dmC);

    // Change in O2 mass [kg]
    const scalar dmO2 = dmC/WC_*Sb_*WO2_;

    // Mass of newly created CO2 [kg]
    const scalar dmCO2 = dmC + dmO2;

    // Update local particle C mass
    dMassSolid[CsLocalId_] += dmC;

    // Update carrier O2 and CO2 mass
    dMassSRCarrier[O2GlobalId_] -= dmO2;
    dMassSRCarrier[CO2GlobalId_] += dmCO2;

    // Heat of reaction [J]
    return -HcCO2_*dmCO2;
}


// ************************************************************************* //
