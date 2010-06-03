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

#include "ReactingMultiphaseParcel.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::GAS(0);

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::LIQ(1);

template<class ParcelType>
const Foam::label Foam::ReactingMultiphaseParcel<ParcelType>::SLD(2);


// * * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::cpEff
(
    TrackData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*td.cloud().composition().cp(idG, YGas_, p, T)
      + this->Y_[LIQ]*td.cloud().composition().cp(idL, YLiquid_, p, T)
      + this->Y_[SLD]*td.cloud().composition().cp(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::HEff
(
    TrackData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*td.cloud().composition().H(idG, YGas_, p, T)
      + this->Y_[LIQ]*td.cloud().composition().H(idL, YLiquid_, p, T)
      + this->Y_[SLD]*td.cloud().composition().H(idS, YSolid_, p, T);
}


template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::LEff
(
    TrackData& td,
    const scalar p,
    const scalar T,
    const label idG,
    const label idL,
    const label idS
) const
{
    return
        this->Y_[GAS]*td.cloud().composition().L(idG, YGas_, p, T)
      + this->Y_[LIQ]*td.cloud().composition().L(idL, YLiquid_, p, T)
      + this->Y_[SLD]*td.cloud().composition().L(idS, YSolid_, p, T);
}


template<class ParcelType>
Foam::scalar Foam::ReactingMultiphaseParcel<ParcelType>::updateMassFractions
(
    const scalar mass0,
    const scalarField& dMassGas,
    const scalarField& dMassLiquid,
    const scalarField& dMassSolid
)
{
    scalarField& YMix = this->Y_;

    scalar massGas =
        this->updateMassFraction(mass0*YMix[GAS], dMassGas, YGas_);
    scalar massLiquid =
        this->updateMassFraction(mass0*YMix[LIQ], dMassLiquid, YLiquid_);
    scalar massSolid =
        this->updateMassFraction(mass0*YMix[SLD], dMassSolid, YSolid_);

    scalar massNew = max(massGas + massLiquid + massSolid, ROOTVSMALL);

    YMix[GAS] = massGas/massNew;
    YMix[LIQ] = massLiquid/massNew;
    YMix[SLD] = 1.0 - YMix[GAS] - YMix[LIQ];

    return massNew;
}


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::ReactingMultiphaseParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ReactingParcel<ParcelType>::setCellValues(td, dt, cellI);
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingMultiphaseParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    scalar massCell = this->massCell(cellI);

    scalar addedMass = 0.0;
    forAll(td.cloud().rhoTrans(), i)
    {
        addedMass += td.cloud().rhoTrans(i)[cellI];
    }

    this->rhoc_ += addedMass/td.cloud().pMesh().cellVolumes()[cellI];

    scalar massCellNew = massCell + addedMass;
    this->Uc_ += td.cloud().UTrans()[cellI]/massCellNew;

    scalar cpEff = 0;
    if (addedMass > ROOTVSMALL)
    {
        forAll(td.cloud().rhoTrans(), i)
        {
            scalar Y = td.cloud().rhoTrans(i)[cellI]/addedMass;
            cpEff +=
                Y*td.cloud().mcCarrierThermo().speciesData()[i].Cp(this->Tc_);
        }
    }
    const scalar cpc = td.cpInterp().psi()[cellI];
    this->cpc_ = (massCell*cpc + addedMass*cpEff)/massCellNew;

    this->Tc_ += td.cloud().hsTrans()[cellI]/(this->cpc_*massCellNew);
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingMultiphaseParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    // Define local properties at beginning of timestep
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;
    const scalar rho0 = this->rho_;
    const scalar T0 = this->T_;
    const scalar cp0 = this->cp_;
    const scalar mass0 = this->mass();

    const scalar pc = this->pc_;

    const scalarField& YMix = this->Y_;
    const label idG = td.cloud().composition().idGas();
    const label idL = td.cloud().composition().idLiquid();
    const label idS = td.cloud().composition().idSolid();


    // Calc surface values
    // ~~~~~~~~~~~~~~~~~~~
    scalar Ts, rhos, mus, Pr, kappa;
    ThermoParcel<ParcelType>::
        calcSurfaceValues(td, cellI, T0, Ts, rhos, mus, Pr, kappa);

    // Reynolds number
    scalar Re = this->Re(U0, d0, rhos, mus);


    // Sources
    //~~~~~~~~

    // Explicit momentum source for particle
    vector Su = vector::zero;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = vector::zero;

    // Explicit enthalpy source for particle
    scalar Sh = 0.0;

    // Sensible enthalpy transfer from the particle to the carrier phase
    scalar dhsTrans = 0.0;


    // Phase change in liquid phase
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Mass transfer due to phase change
    scalarField dMassPC(YLiquid_.size(), 0.0);

    // Molar flux of species emitted from the particle (kmol/m^2/s)
    scalar Ne = 0.0;

    // Sum Ni*Cpi*Wi of emission species
    scalar NCpW = 0.0;

    // Surface concentrations of emitted species
    scalarField Cs(td.cloud().mcCarrierThermo().species().size(), 0.0);

    // Calc mass and enthalpy transfer due to phase change
    calcPhaseChange
    (
        td,
        dt,
        cellI,
        Re,
        Ts,
        mus/rhos,
        d0,
        T0,
        mass0,
        idL,
        YMix[LIQ],
        YLiquid_,
        dMassPC,
        Sh,
        Ne,
        NCpW,
        Cs
    );


    // Devolatilisation
    // ~~~~~~~~~~~~~~~~

    // Mass transfer due to devolatilisation
    scalarField dMassDV(YGas_.size(), 0.0);

    // Calc mass and enthalpy transfer due to devolatilisation
    calcDevolatilisation
    (
        td,
        dt,
        Ts,
        d0,
        T0,
        mass0,
        this->mass0_,
        idG,
        YMix[GAS],
        YGas_,
        canCombust_,
        dMassDV,
        Sh,
        Ne,
        NCpW,
        Cs
    );

    // Correct surface values due to emitted species
    correctSurfaceValues(td, cellI, Ts, Cs, rhos, mus, Pr, kappa);


    // Surface reactions
    // ~~~~~~~~~~~~~~~~~

    // Change in carrier phase composition due to surface reactions
    scalarField dMassSRGas(YGas_.size(), 0.0);
    scalarField dMassSRLiquid(YLiquid_.size(), 0.0);
    scalarField dMassSRSolid(YSolid_.size(), 0.0);
    scalarField
        dMassSRCarrier
        (
            td.cloud().mcCarrierThermo().species().size(),
            0.0
        );

    // Clac mass and enthalpy transfer due to surface reactions
    calcSurfaceReactions
    (
        td,
        dt,
        cellI,
        d0,
        T0,
        mass0,
        canCombust_,
        Ne,
        YMix,
        YGas_,
        YLiquid_,
        YSolid_,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier,
        Sh,
        dhsTrans
    );


    // Update component mass fractions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalarField dMassGas = dMassDV + dMassSRGas;
    scalarField dMassLiquid = dMassPC + dMassSRLiquid;
    scalarField dMassSolid = dMassSRSolid;

    scalar mass1 =
        updateMassFractions(mass0, dMassGas, dMassLiquid, dMassSolid);


    // Heat transfer
    // ~~~~~~~~~~~~~

    // Calculate new particle temperature
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

    dUTrans += 0.5*(mass0 - mass1)*(U0 + U1);


    // Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (td.cloud().coupled())
    {
        // Transfer mass lost from particle to carrier mass source
        forAll(YGas_, i)
        {
            label gid = td.cloud().composition().localToGlobalCarrierId(GAS, i);
            td.cloud().rhoTrans(gid)[cellI] += np0*dMassGas[i];
        }
        forAll(YLiquid_, i)
        {
            label gid = td.cloud().composition().localToGlobalCarrierId(LIQ, i);
            td.cloud().rhoTrans(gid)[cellI] += np0*dMassLiquid[i];
        }
/*
        // No mapping between solid components and carrier phase
        forAll(YSolid_, i)
        {
            label gid = td.cloud().composition().localToGlobalCarrierId(SLD, i);
            td.cloud().rhoTrans(gid)[cellI] += np0*dMassSolid[i];
        }
*/
        forAll(dMassSRCarrier, i)
        {
            td.cloud().rhoTrans(i)[cellI] += np0*dMassSRCarrier[i];
        }

        // Update momentum transfer
        td.cloud().UTrans()[cellI] += np0*dUTrans;

        // Update sensible enthalpy transfer
        td.cloud().hsTrans()[cellI] += np0*dhsTrans;
    }


    // Remove the particle when mass falls below minimum threshold
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (mass1 < td.constProps().minParticleMass())
    {
        td.keepParticle = false;

        if (td.cloud().coupled())
        {
            // Absorb parcel into carrier phase
            forAll(YGas_, i)
            {
                label gid =
                    td.cloud().composition().localToGlobalCarrierId(GAS, i);
                td.cloud().rhoTrans(gid)[cellI] += np0*mass1*YMix[GAS]*YGas_[i];
            }
            forAll(YLiquid_, i)
            {
                label gid =
                    td.cloud().composition().localToGlobalCarrierId(LIQ, i);
                td.cloud().rhoTrans(gid)[cellI] +=
                    np0*mass1*YMix[LIQ]*YLiquid_[i];
            }
/*
            // No mapping between solid components and carrier phase
            forAll(YSolid_, i)
            {
                label gid =
                    td.cloud().composition().localToGlobalCarrierId(SLD, i);
                td.cloud().rhoTrans(gid)[cellI] +=
                    np0*mass1*YMix[SLD]*YSolid_[i];
            }
*/
            td.cloud().UTrans()[cellI] += np0*mass1*U1;
            td.cloud().hsTrans()[cellI] +=
                np0*mass1*HEff(td, pc, T1, idG, idL, idS); // using total h
        }
    }


    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    else
    {
        this->cp_ = cpEff(td, pc, T1, idG, idL, idS);
        this->T_ = T1;
        this->U_ = U1;

        // Update particle density or diameter
        if (td.constProps().constantVolume())
        {
            this->rho_ = mass1/this->volume();
        }
        else
        {
            this->d_ = cbrt(mass1/this->rho_*6.0/mathematicalConstant::pi);
        }
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcDevolatilisation
(
    TrackData& td,
    const scalar dt,
    const scalar Ts,
    const scalar d,
    const scalar T,
    const scalar mass,
    const scalar mass0,
    const label idVolatile,
    const scalar YVolatileTot,
    const scalarField& YVolatile,
    bool& canCombust,
    scalarField& dMassDV,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
) const
{
    // Check that model is active, and that the parcel temperature is
    // within necessary limits for devolatilisation to occur
    if
    (
        !td.cloud().devolatilisation().active()
     || T < td.constProps().Tvap()
    )
    {
        return;
    }

    // Total mass of volatiles evolved
    const scalar dMassTot = td.cloud().devolatilisation().calculate
    (
        dt,
        mass0,
        mass,
        T,
        td.cloud().composition().YMixture0()[idVolatile],
        YVolatileTot,
        canCombust
    );

    // Volatile mass transfer - equal components of each volatile specie
    dMassDV = YVolatile*dMassTot;

    td.cloud().addToMassDevolatilisation(this->nParticle_*dMassTot);

    Sh -= dMassTot*td.constProps().LDevol()/dt;

    // Molar average molecular weight of carrier mix
    const scalar Wc = this->rhoc_*specie::RR*this->Tc_/this->pc_;

    // Update molar emissions
    forAll(dMassDV, i)
    {
        // Note: hardcoded gaseous diffusivities for now
        // TODO: add to carrier thermo
        const scalar beta = sqr(cbrt(15.0) + cbrt(15.0));
        const label id =
            td.cloud().composition().localToGlobalCarrierId(GAS, i);
        const scalar Cp = td.cloud().mcCarrierThermo().speciesData()[id].Cp(Ts);
        const scalar W = td.cloud().mcCarrierThermo().speciesData()[id].W();
        const scalar Ni = dMassDV[i]/(this->areaS(d)*dt*W);

        // Dab calc'd using API vapour mass diffusivity function
        const scalar Dab =
            3.6059e-3*(pow(1.8*Ts, 1.75))*sqrt(1.0/W + 1.0/Wc)/(this->pc_*beta);

        N += Ni;
        NCpW += Ni*Cp*W;
        Cs[id] += Ni*d/(2.0*Dab);
     }
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingMultiphaseParcel<ParcelType>::calcSurfaceReactions
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar d,
    const scalar T,
    const scalar mass,
    const bool canCombust,
    const scalar N,
    const scalarField& YMix,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    scalarField& dMassSRGas,
    scalarField& dMassSRLiquid,
    scalarField& dMassSRSolid,
    scalarField& dMassSRCarrier,
    scalar& Sh,
    scalar& dhsTrans
) const
{
    // Check that model is active
    if (!td.cloud().surfaceReaction().active() || !canCombust)
    {
        return;
    }

    // Update surface reactions
    const scalar hReaction = td.cloud().surfaceReaction().calculate
    (
        dt,
        cellI,
        d,
        T,
        this->Tc_,
        this->pc_,
        this->rhoc_,
        mass,
        YGas,
        YLiquid,
        YSolid,
        YMix,
        N,
        dMassSRGas,
        dMassSRLiquid,
        dMassSRSolid,
        dMassSRCarrier
    );

    td.cloud().addToMassSurfaceReaction
    (
        this->nParticle_
       *(sum(dMassSRGas) + sum(dMassSRLiquid) + sum(dMassSRSolid))
    );

    const scalar xsi = min(T/5000.0, 1.0);
    const scalar hRetentionCoeffMod =
        (1.0 - xsi*xsi)*td.constProps().hRetentionCoeff();

    Sh += hRetentionCoeffMod*hReaction/dt;

    dhsTrans += (1.0 - hRetentionCoeffMod)*hReaction;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ParcelType>
Foam::ReactingMultiphaseParcel<ParcelType>::ReactingMultiphaseParcel
(
    const ReactingMultiphaseParcel<ParcelType>& p
)
:
    ReactingParcel<ParcelType>(p),
    YGas_(p.YGas_),
    YLiquid_(p.YLiquid_),
    YSolid_(p.YSolid_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingMultiphaseParcelIO.C"

// ************************************************************************* //
