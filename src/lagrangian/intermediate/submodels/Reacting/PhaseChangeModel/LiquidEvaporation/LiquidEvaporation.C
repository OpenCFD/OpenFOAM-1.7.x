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

#include "LiquidEvaporation.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalarField Foam::LiquidEvaporation<CloudType>::calcXc
(
    const label cellI
) const
{
    scalarField Xc(this->owner().mcCarrierThermo().Y().size());

    scalar Winv = 0.0;
    forAll(Xc, i)
    {
        scalar Y = this->owner().mcCarrierThermo().Y()[i][cellI];
        Winv += Y/this->owner().mcCarrierThermo().speciesData()[i].W();
        Xc[i] = Y/this->owner().mcCarrierThermo().speciesData()[i].W();
    }

    return Xc/Winv;
}


template <class CloudType>
Foam::scalar Foam::LiquidEvaporation<CloudType>::Sh
(
    const scalar Re,
    const scalar Sc
) const
{
    return 2.0 + 0.6*Foam::sqrt(Re)*pow(Sc, 0.333);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::LiquidEvaporation<CloudType>::LiquidEvaporation
(
    const dictionary& dict,
    CloudType& owner
)
:
    PhaseChangeModel<CloudType>(dict, owner, typeName),
    liquids_
    (
        liquidMixture::New
        (
            owner.mesh().objectRegistry::lookupObject<dictionary>
            (
                owner.carrierThermo().name()
            )
        )
    ),
    activeLiquids_(this->coeffDict().lookup("activeLiquids")),
    liqToCarrierMap_(activeLiquids_.size(), -1),
    liqToLiqMap_(activeLiquids_.size(), -1)
{
    if (activeLiquids_.size() == 0)
    {
        WarningIn
        (
            "Foam::LiquidEvaporation<CloudType>::LiquidEvaporation"
            "("
                "const dictionary& dict, "
                "CloudType& owner"
            ")"
        )   << "Evaporation model selected, but no active liquids defined"
            << nl << endl;
    }

    // Determine mapping between liquid and carrier phase species
    forAll(activeLiquids_, i)
    {
        liqToCarrierMap_[i] =
            owner.composition().globalCarrierId(activeLiquids_[i]);
    }

    // Determine mapping between model active liquids and global liquids
    label idLiquid = owner.composition().idLiquid();
    forAll(activeLiquids_, i)
    {
        liqToLiqMap_[i] =
            owner.composition().localId(idLiquid, activeLiquids_[i]);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::LiquidEvaporation<CloudType>::~LiquidEvaporation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::LiquidEvaporation<CloudType>::active() const
{
    return true;
}


template<class CloudType>
void Foam::LiquidEvaporation<CloudType>::calculate
(
    const scalar dt,
    const label cellI,
    const scalar d,
    const scalar T,
    const scalar pc,
    const scalar Tc,
    const scalar nuc,
    const vector& Ur,
    scalarField& dMassPC
) const
{
    // construct carrier phase species volume fractions for cell, cellI
    scalarField Xc = calcXc(cellI);

    // droplet surface area
    scalar A = mathematicalConstant::pi*sqr(d);

    // Reynolds number
    scalar Re = mag(Ur)*d/(nuc + ROOTVSMALL);

    // film temperature evaluated using the 2/3 rule
    scalar Tf = (2.0*T + Tc)/3.0;

    // calculate mass transfer of each specie in liquid
    forAll(activeLiquids_, i)
    {
        label gid = liqToCarrierMap_[i];
        label lid = liqToLiqMap_[i];

        // vapour diffusivity at film temperature and cell pressure [m2/s]
        scalar Dab = liquids_->properties()[lid].D(pc, Tf);

        // saturation pressure for species i at film temperature and cell
        // pressure [pa] - carrier phase pressure assumed equal to the liquid
        // vapour pressure close to the surface
        // - limited to pc if pSat > pc
        scalar pSat = min(liquids_->properties()[lid].pv(pc, Tf), pc);

        // Schmidt number
        scalar Sc = nuc/(Dab + ROOTVSMALL);

        // Sherwood number
        scalar Sh = this->Sh(Re, Sc);

        // mass transfer coefficient [m/s]
        scalar kc = Sh*Dab/(d + ROOTVSMALL);

        // vapour concentration at droplet surface at film temperature [kmol/m3]
        scalar Cs = pSat/(specie::RR*Tf);

        // vapour concentration in bulk gas at film temperature [kmol/m3]
        scalar Cinf = Xc[gid]*pc/(specie::RR*Tf);

        // molar flux of vapour [kmol/m2/s]
        scalar Ni = max(kc*(Cs - Cinf), 0.0);

        // mass transfer [kg]
        dMassPC[lid] += Ni*A*liquids_->properties()[lid].W()*dt;
    }
}


// ************************************************************************* //
