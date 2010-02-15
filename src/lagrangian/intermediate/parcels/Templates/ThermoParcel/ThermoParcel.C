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

#include "ThermoParcel.H"
#include "radiationConstants.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::ThermoParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    KinematicParcel<ParcelType>::setCellValues(td, dt, cellI);

    cpc_ = td.cpInterp().interpolate(this->position(), cellI);

    Tc_ = td.TInterp().interpolate(this->position(), cellI);

    if (Tc_ < td.constProps().TMin())
    {
        WarningIn
        (
            "void Foam::ThermoParcel<ParcelType>::setCellValues"
            "("
                "TrackData&, "
                "const scalar, "
                "const label"
            ")"
        )   << "Limiting observed temperature in cell " << cellI << " to "
            << td.constProps().TMin() <<  nl << endl;

        Tc_ = td.constProps().TMin();
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ThermoParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    this->Uc_ += td.cloud().UTrans()[cellI]/this->massCell(cellI);

    scalar cpMean = td.cpInterp().psi()[cellI];
    Tc_ += td.cloud().hsTrans()[cellI]/(cpMean*this->massCell(cellI));
}


template<class ParcelType>
template<class TrackData>
void Foam::ThermoParcel<ParcelType>::calcSurfaceValues
(
    TrackData& td,
    const label cellI,
    const scalar T,
    scalar& Ts,
    scalar& rhos,
    scalar& mus,
    scalar& Pr,
    scalar& kappa
) const
{
    // Surface temperature using two thirds rule
    Ts = (2.0*T + Tc_)/3.0;

    // Assuming thermo props vary linearly with T for small dT
    scalar factor = td.TInterp().interpolate(this->position(), cellI)/Ts;
    rhos = this->rhoc_*factor;
    mus = td.muInterp().interpolate(this->position(), cellI)/factor;

    Pr = td.constProps().Pr();
    kappa = cpc_*mus/Pr;
}


template<class ParcelType>
template<class TrackData>
void Foam::ThermoParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector U0 = this->U_;
    const scalar rho0 = this->rho_;
    const scalar T0 = this->T_;
    const scalar cp0 = this->cp_;
    const scalar mass0 = this->mass();


    // Calc surface values
    // ~~~~~~~~~~~~~~~~~~~
    scalar Ts, rhos, mus, Pr, kappa;
    calcSurfaceValues(td, cellI, T0, Ts, rhos, mus, Pr, kappa);

    // Reynolds number
    scalar Re = this->Re(U0, d0, rhos, mus);


    // Sources
    // ~~~~~~~

    // Explicit momentum source for particle
    vector Su = vector::zero;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;

    // Explicit enthalpy source for particle
    scalar Sh = 0.0;

    // Sensible enthalpy transfer from the particle to the carrier phase
    scalar dhsTrans = 0.0;


    // Heat transfer
    // ~~~~~~~~~~~~~

    // Sum Ni*Cpi*Wi of emission species
    scalar NCpW = 0.0;

    // Calculate new particle velocity
    scalar T1 =
        calcHeatTransfer
        (
            td,
            dt,
            cellI,
            Re,
            Pr,
            kappa,
            d0,
            rho0,
            T0,
            cp0,
            NCpW,
            Sh,
            dhsTrans
        );


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    vector U1 =
        calcVelocity(td, dt, cellI, Re, mus, d0, U0, rho0, mass0, Su, dUTrans);


    //  Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (td.cloud().coupled())
    {
        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*dUTrans;

        // Update sensible enthalpy transfer
        td.cloud().hsTrans()[cellI] += np0*dhsTrans;
    }

    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    this->U_ = U1;
    T_ = T1;
}


template<class ParcelType>
template <class TrackData>
Foam::scalar Foam::ThermoParcel<ParcelType>::calcHeatTransfer
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar Re,
    const scalar Pr,
    const scalar kappa,
    const scalar d,
    const scalar rho,
    const scalar T,
    const scalar cp,
    const scalar NCpW,
    const scalar Sh,
    scalar& dhsTrans
)
{
    if (!td.cloud().heatTransfer().active())
    {
        return T;
    }

    // Calc heat transfer coefficient
    scalar htc = td.cloud().heatTransfer().htc(d, Re, Pr, kappa, NCpW);

    if (mag(htc) < ROOTVSMALL && !td.cloud().radiation())
    {
        return max(T + dt*Sh/(this->volume(d)*rho*cp), td.constProps().TMin());
    }

    const scalar As = this->areaS(d);
    scalar ap = Tc_ + Sh/As/htc;
    scalar bp = 6.0*(Sh/As + htc*(Tc_ - T));
    if (td.cloud().radiation())
    {
        const scalarField& G =
            td.cloud().mesh().objectRegistry::lookupObject<volScalarField>("G");
        const scalar Gc = G[cellI];
        const scalar sigma = radiation::sigmaSB.value();
        const scalar epsilon = td.constProps().epsilon0();

        ap = (ap + epsilon*Gc/(4.0*htc))/(1.0 + epsilon*sigma*pow3(T)/htc);
        bp += 6.0*(epsilon*(Gc/4.0 - sigma*pow4(T)));
    }
    bp /= rho*d*cp*(ap - T);

    // Integrate to find the new parcel temperature
    IntegrationScheme<scalar>::integrationResult Tres =
        td.cloud().TIntegrator().integrate(T, dt, ap, bp);

    scalar Tnew = max(Tres.value(), td.constProps().TMin());

    dhsTrans += dt*htc*As*(0.5*(T + Tnew) - Tc_);

    return Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ParcelType>
Foam::ThermoParcel<ParcelType>::ThermoParcel
(
    const ThermoParcel<ParcelType>& p
)
:
    KinematicParcel<ParcelType>(p),
    T_(p.T_),
    cp_(p.cp_),
    Tc_(p.Tc_),
    cpc_(p.cpc_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ThermoParcelIO.C"

// ************************************************************************* //
