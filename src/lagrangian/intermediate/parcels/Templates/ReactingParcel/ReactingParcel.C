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

    // Calc mass and enthalpy transfer due to phase change
    calcPhaseChange
    (
        td,
        dt,
        cellI,
        d0,
        T0,
        U0,
        mass0,
        0,
        1.0,
        Y_,
        dMassPC,
        Sh,
        dhsTrans
    );

    // Update particle component mass and mass fractions
    scalar mass1 = updateMassFraction(mass0, dMassPC, Y_);


    // Heat transfer
    // ~~~~~~~~~~~~~

    // Calculate new particle temperature
    scalar T1 =
        calcHeatTransfer(td, dt, cellI, d0, U0, rho0, T0, cp0, Sh, dhsTrans);


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    vector U1 = calcVelocity(td, dt, cellI, d0, U0, rho0, mass0, Su, dUTrans);

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
            td.cloud().hcTrans()[cellI] +=
                np0
               *dMassPC[i]
               *td.cloud().mcCarrierThermo().speciesData()[gid].H(T0);
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
            td.cloud().hcTrans()[cellI] +=
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
    const scalar d,
    const scalar T,
    const vector& U,
    const scalar mass,
    const label idPhase,
    const scalar YPhase,
    const scalarField& YComponents,
    scalarField& dMassPC,
    scalar& Sh,
    scalar& dhsTrans
)
{
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
        d,
        min(T, td.constProps().Tbp()), // Limit to boiling temperature
        pc_,
        this->Tc_,
        this->muc_/(this->rhoc_ + ROOTVSMALL),
        U - this->Uc_,
        dMassPC
    );

    // Limit phase change mass by availability of each specie
    dMassPC = min(mass*YPhase*YComponents, dMassPC);

    scalar dMassTot = sum(dMassPC);

    // Add to cumulative phase change mass
    td.cloud().addToMassPhaseChange(this->nParticle_*dMassTot);

    // Enthalphy transfer to carrier phase
    label id;
    forAll(YComponents, i)
    {
        id = td.cloud().composition().localToGlobalCarrierId(idPhase, i);
        const scalar hv = td.cloud().mcCarrierThermo().speciesData()[id].H(T);

        id = td.cloud().composition().globalIds(idPhase)[i];
        const scalar hl =
            td.cloud().composition().liquids().properties()[id].h(pc_, T);

        Sh += dMassPC[i]*(hl - hv)/dt;
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

