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

#include "COxidationMurphyShaddix.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::COxidationMurphyShaddix<CloudType>::maxIters_ = 1000;

template<class CloudType>
Foam::scalar Foam::COxidationMurphyShaddix<CloudType>::tolerance_ = 1e-06;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationMurphyShaddix<CloudType>::COxidationMurphyShaddix
(
    const dictionary& dict,
    CloudType& owner
)
:
    SurfaceReactionModel<CloudType>
    (
        dict,
        owner,
        typeName
    ),
    D0_(dimensionedScalar(this->coeffDict().lookup("D0")).value()),
    rho0_(dimensionedScalar(this->coeffDict().lookup("rho0")).value()),
    T0_(dimensionedScalar(this->coeffDict().lookup("T0")).value()),
    Dn_(dimensionedScalar(this->coeffDict().lookup("Dn")).value()),
    A_(dimensionedScalar(this->coeffDict().lookup("A")).value()),
    E_(dimensionedScalar(this->coeffDict().lookup("E")).value()),
    n_(dimensionedScalar(this->coeffDict().lookup("n")).value()),
    WVol_(dimensionedScalar(this->coeffDict().lookup("WVol")).value()),
    CsLocalId_(-1),
    O2GlobalId_(owner.composition().globalCarrierId("O2")),
    CO2GlobalId_(owner.composition().globalCarrierId("CO2")),
    WC_(0.0),
    WO2_(0.0)
{
    // Determine Cs ids
    label idSolid = owner.composition().idSolid();
    CsLocalId_ = owner.composition().localId(idSolid, "C");

    // Set local copies of thermo properties
    WO2_ = owner.mcCarrierThermo().speciesData()[O2GlobalId_].W();
    scalar WCO2 = owner.mcCarrierThermo().speciesData()[CO2GlobalId_].W();
    WC_ = WCO2 - WO2_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::COxidationMurphyShaddix<CloudType>::~COxidationMurphyShaddix()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::COxidationMurphyShaddix<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::COxidationMurphyShaddix<CloudType>::calculate
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
    const scalar N,
    scalarField& dMassGas,
    scalarField& dMassLiquid,
    scalarField& dMassSolid,
    scalarField& dMassSRCarrier
) const
{
    // Fraction of remaining combustible material
    const label idSolid = CloudType::parcelType::SLD;
    const scalar fComb = YMixture[idSolid]*YSolid[CsLocalId_];

    // Surface combustion until combustible fraction is consumed
    if (fComb < SMALL)
    {
        return 0.0;
    }

    // Cell carrier phase O2 species density [kg/m^3]
    const scalar rhoO2 =
        rhoc*this->owner().mcCarrierThermo().Y(O2GlobalId_)[cellI];

    if (rhoO2 < SMALL)
    {
        return 0.0;
    }

    // Particle surface area [m^2]
    const scalar Ap = mathematicalConstant::pi*sqr(d);

    // Calculate diffision constant at continuous phase temperature
    // and density [m^2/s]
    const scalar D = D0_*(rho0_/rhoc)*pow(Tc/T0_, Dn_);

    // Far field partial pressure O2 [Pa]
    const scalar ppO2 = rhoO2/WO2_*specie::RR*Tc;

    // Total molar concentration of the carrier phase [kmol/m^3]
    const scalar C = pc/(specie::RR*Tc);

    if (debug)
    {
        Pout<< "mass  = " << mass << nl
            << "fComb = " << fComb << nl
            << "Ap    = " << Ap << nl
            << "dt    = " << dt << nl
            << "C     = " << C << nl
            << endl;
    }

    // Molar reaction rate per unit surface area [kmol/(m^2.s)]
    scalar qCsOld = 0;
    scalar qCs = 1;

    const scalar qCsLim = mass*fComb/(WC_*Ap*dt);

    if (debug)
    {
        Pout << "qCsLim = " << qCsLim << endl;
    }

    label iter = 0;
    while ((mag(qCs - qCsOld)/qCs > tolerance_) && (iter <= maxIters_))
    {
        qCsOld = qCs;
        const scalar PO2Surface = ppO2*exp(-(qCs + N)*d/(2*C*D));
        qCs = A_*exp(-E_/(specie::RR*T))*pow(PO2Surface, n_);
        qCs = (100.0*qCs + iter*qCsOld)/(100.0 + iter);
        qCs = min(qCs, qCsLim);

        if (debug)
        {
            Pout<< "iter = " << iter
                << ", qCsOld = " << qCsOld
                << ", qCs = " << qCs
                << nl << endl;
        }

        iter++;
    }

    if (iter > maxIters_)
    {
        WarningIn
        (
            "scalar Foam::COxidationMurphyShaddix<CloudType>::calculate(...)"
        )   << "iter limit reached (" << maxIters_ << ")" << nl << endl;
    }

    // Calculate the number of molar units reacted
    scalar dOmega = qCs*Ap*dt;

    // Add to carrier phase mass transfer
    dMassSRCarrier[O2GlobalId_] += -dOmega*WO2_;
    dMassSRCarrier[CO2GlobalId_] += dOmega*(WC_ + WO2_);

    // Add to particle mass transfer
    dMassSolid[CsLocalId_] += dOmega*WC_;

    const scalar HC =
        this->owner().composition().solids().properties()[CsLocalId_].Hf()
      + this->owner().composition().solids().properties()[CsLocalId_].cp()*T;
    const scalar HCO2 =
        this->owner().mcCarrierThermo().speciesData()[CO2GlobalId_].H(T);
    const scalar HO2 =
        this->owner().mcCarrierThermo().speciesData()[O2GlobalId_].H(T);

    // Heat of reaction
    return dOmega*(WC_*HC + WO2_*HO2 - (WC_ + WO2_)*HCO2);
}


// ************************************************************************* //
