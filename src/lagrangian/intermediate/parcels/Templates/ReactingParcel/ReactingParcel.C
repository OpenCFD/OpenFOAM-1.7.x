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

#include "ReactingParcel.H"
#include "mathematicalConstants.H"
#include "specie.H"

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ThermoParcel<ParcelType>::setCellValues(td, dt, cellI);

    pc_ = td.pInterp().interpolate(this->position(), cellI);
    if (pc_ < td.constProps().pMin())
    {
        WarningIn
        (
            "void Foam::ReactingParcel<ParcelType>::setCellValues"
            "("
                "TrackData&, "
                "const scalar, "
                "const label"
            ")"
        )   << "Limiting observed pressure in cell " << cellI << " to "
            << td.constProps().pMin() <<  nl << endl;

        pc_ = td.constProps().pMin();
    }
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::cellValueSourceCorrection
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
void Foam::ReactingParcel<ParcelType>::correctSurfaceValues
(
    TrackData& td,
    const label cellI,
    const scalar T,
    const scalarField& Cs,
    scalar& rhos,
    scalar& mus,
    scalar& Pr,
    scalar& kappa
)
{
    // No correction if total concentration of emitted species is small
    if (sum(Cs) < SMALL)
    {
        return;
    }

    // Far field carrier  molar fractions
    scalarField Xinf(td.cloud().mcCarrierThermo().speciesData().size());

    forAll(Xinf, i)
    {
        Xinf[i] =
            td.cloud().mcCarrierThermo().Y(i)[cellI]
           /td.cloud().mcCarrierThermo().speciesData()[i].W();
    }
    Xinf /= sum(Xinf);

    // Molar fraction of far field species at particle surface
    const scalar Xsff = 1.0 - min(sum(Cs)*specie::RR*this->T_/pc_, 1.0);

    // Surface carrier total molar concentration
    const scalar CsTot = pc_/(specie::RR*this->T_);

    // Surface carrier composition (molar fraction)
    scalarField Xs(Xinf.size());

    // Surface carrier composition (mass fraction)
    scalarField Ys(Xinf.size());

    forAll(Xs, i)
    {
        // Molar concentration of species at particle surface
        const scalar Csi = Cs[i] + Xsff*Xinf[i]*CsTot;

        Xs[i] = (2.0*Csi + Xinf[i]*CsTot)/3.0;
        Ys[i] = Xs[i]*td.cloud().mcCarrierThermo().speciesData()[i].W();
    }
    Xs /= sum(Xs);
    Ys /= sum(Ys);


    rhos = 0;
    mus = 0;
    kappa = 0;
    scalar cps = 0;
    scalar sumYiSqrtW = 0;
    scalar sumYiCbrtW = 0;

    forAll(Ys, i)
    {
        const scalar sqrtW =
            sqrt(td.cloud().mcCarrierThermo().speciesData()[i].W());
        const scalar cbrtW =
            cbrt(td.cloud().mcCarrierThermo().speciesData()[i].W());

        rhos += Xs[i]*td.cloud().mcCarrierThermo().speciesData()[i].W();
        mus += Ys[i]*sqrtW*td.cloud().mcCarrierThermo().speciesData()[i].mu(T);
        kappa +=
            Ys[i]*cbrtW*td.cloud().mcCarrierThermo().speciesData()[i].kappa(T);
        cps += Xs[i]*td.cloud().mcCarrierThermo().speciesData()[i].Cp(T);

        sumYiSqrtW += Ys[i]*sqrtW;
        sumYiCbrtW += Ys[i]*cbrtW;
    }

    rhos *= pc_/(specie::RR*T);
    mus /= sumYiSqrtW;
    kappa /= sumYiCbrtW;
    Pr = cps*mus/kappa;
}


template<class ParcelType>
Foam::scalar Foam::ReactingParcel<ParcelType>::updateMassFraction
(
    const scalar mass0,
    const scalarField& dMass,
    scalarField& Y
) const
{
    scalar mass1 = mass0 - sum(dMass);

    // only update the mass fractions if the new particle mass is finite
    if (mass1 > ROOTVSMALL)
    {
        forAll(Y, i)
        {
            Y[i] = (Y[i]*mass0 - dMass[i])/mass1;
        }
    }

    return mass1;
}


template<class ParcelType>
template<class TrackData>
void Foam::ReactingParcel<ParcelType>::calc
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
    const vector& U0 = this->U_;
    const scalar rho0 = this->rho_;
    const scalar T0 = this->T_;
    const scalar cp0 = this->cp_;
    const scalar mass0 = this->mass();


    // Calc surface values
    // ~~~~~~~~~~~~~~~~~~~
    scalar Ts, rhos, mus, Pr, kappa;
    this->calcSurfaceValues(td, cellI, T0, Ts, rhos, mus, Pr, kappa);

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


    // Phase change
    // ~~~~~~~~~~~~

    // Mass transfer due to phase change
    scalarField dMassPC(Y_.size(), 0.0);

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
        0,
        1.0,
        Y_,
        dMassPC,
        Sh,
        Ne,
        NCpW,
        Cs
    );

    // Correct surface values due to emitted species
    correctSurfaceValues(td, cellI, Ts, Cs, rhos, mus, Pr, kappa);

    // Update particle component mass and mass fractions
    scalar mass1 = updateMassFraction(mass0, dMassPC, Y_);


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
        forAll(dMassPC, i)
        {
            label gid = td.cloud().composition().localToGlobalCarrierId(0, i);
            td.cloud().rhoTrans(gid)[cellI] += np0*dMassPC[i];
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
            forAll(Y_, i)
            {
                label gid =
                    td.cloud().composition().localToGlobalCarrierId(0, i);
                td.cloud().rhoTrans(gid)[cellI] += np0*mass1*Y_[i];
            }
            td.cloud().UTrans()[cellI] += np0*mass1*U1;
            td.cloud().hsTrans()[cellI] +=
                np0*mass1*td.cloud().composition().H(0, Y_, pc_, T1);
        }
    }


    // Set new particle properties
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    else
    {
        this->cp_ = td.cloud().composition().cp(0, Y_, pc_, T1);
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
void Foam::ReactingParcel<ParcelType>::calcPhaseChange
(
    TrackData& td,
    const scalar dt,
    const label cellI,
    const scalar Re,
    const scalar Ts,
    const scalar nus,
    const scalar d,
    const scalar T,
    const scalar mass,
    const label idPhase,
    const scalar YPhase,
    const scalarField& YComponents,
    scalarField& dMassPC,
    scalar& Sh,
    scalar& N,
    scalar& NCpW,
    scalarField& Cs
)
{
    typedef PhaseChangeModel
    <
        typename ReactingParcel<ParcelType>::trackData::cloudType
    > phaseChangeModelType;

    if
    (
        !td.cloud().phaseChange().active()
     || T < td.constProps().Tvap()
     || YPhase < SMALL
    )
    {
        return;
    }

    // Calculate mass transfer due to phase change
    td.cloud().phaseChange().calculate
    (
        dt,
        cellI,
        Re,
        d,
        nus,
        T,
        Ts,
        pc_,
        dMassPC
    );

    // Limit phase change mass by availability of each specie
    dMassPC = min(mass*YPhase*YComponents, dMassPC);

    scalar dMassTot = sum(dMassPC);

    // Add to cumulative phase change mass
    td.cloud().addToMassPhaseChange(this->nParticle_*dMassTot);

    // Average molecular weight of carrier mix - assumes perfect gas
    scalar Wc = this->rhoc_*specie::RR*this->Tc_/this->pc_;

    forAll(YComponents, i)
    {
        const label idc =
            td.cloud().composition().localToGlobalCarrierId(idPhase, i);
        const label idl = td.cloud().composition().globalIds(idPhase)[i];

        // Calculate enthalpy transfer
        if
        (
            td.cloud().phaseChange().enthalpyTransfer()
         == phaseChangeModelType::etLatentHeat
        )
        {
            scalar hlp =
                td.cloud().composition().liquids().properties()[idl].hl(pc_, T);

            Sh -= dMassPC[i]*hlp/dt;
        }
        else
        {
            // Note: enthalpies of both phases must use the same reference
            scalar hc = td.cloud().mcCarrierThermo().speciesData()[idc].H(T);
            scalar hp =
                td.cloud().composition().liquids().properties()[idl].h(pc_, T);

            Sh -= dMassPC[i]*(hc - hp)/dt;
        }

        // Update particle surface thermo properties
        const scalar Dab =
            td.cloud().composition().liquids().properties()[idl].D(pc_, Ts, Wc);

        const scalar Cp =
            td.cloud().mcCarrierThermo().speciesData()[idc].Cp(Ts);
        const scalar W = td.cloud().mcCarrierThermo().speciesData()[idc].W();
        const scalar Ni = dMassPC[i]/(this->areaS(d)*dt*W);

        // Molar flux of species coming from the particle (kmol/m^2/s)
        N += Ni;

        // Sum of Ni*Cpi*Wi of emission species
        NCpW += Ni*Cp*W;

        // Concentrations of emission species
        Cs[idc] += Ni*d/(2.0*Dab);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ParcelType>
Foam::ReactingParcel<ParcelType>::ReactingParcel
(
    const ReactingParcel<ParcelType>& p
)
:
    ThermoParcel<ParcelType>(p),
    mass0_(p.mass0_),
    Y_(p.Y_),
    pc_(p.pc_)
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ReactingParcelIO.C"

// ************************************************************************* //

