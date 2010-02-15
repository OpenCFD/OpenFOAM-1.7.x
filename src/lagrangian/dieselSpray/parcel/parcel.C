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

#include "parcel.H"

#include "spray.H"
#include "dragModel.H"
#include "evaporationModel.H"
#include "heatTransferModel.H"
#include "wallModel.H"
#include "wallPolyPatch.H"
#include "wedgePolyPatch.H"
#include "processorPolyPatch.H"
#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineParticleTypeNameAndDebug(parcel, 0);
    defineTemplateTypeNameAndDebug(Cloud<parcel>, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parcel::parcel
(
    const Cloud<parcel>& cloud,
    const vector& position,
    const label cellI,
    const vector& n,
    const scalar d,
    const scalar T,
    const scalar m,
    const scalar y,
    const scalar yDot,
    const scalar ct,
    const scalar ms,
    const scalar tTurb,
    const scalar liquidCore,
    const scalar injector,
    const vector& U,
    const vector& Uturb,
    const scalarField& X,
    const List<word>& liquidNames
)
:
    Particle<parcel>(cloud, position, cellI),
    liquidComponents_
    (
        liquidNames
    ),
    d_(d),
    T_(T),
    m_(m),
    y_(y),
    yDot_(yDot),
    ct_(ct),
    ms_(ms),
    tTurb_(tTurb),
    liquidCore_(liquidCore),
    injector_(injector),
    U_(U),
    Uturb_(Uturb),
    n_(n),
    X_(X),
    tMom_(GREAT)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::parcel::move(spray& sDB)
{
    const polyMesh& mesh = cloud().pMesh();
    const polyBoundaryMesh& pbMesh = mesh.boundaryMesh();

    const liquidMixture& fuels = sDB.fuels();

    scalar deltaT = sDB.runTime().deltaT().value();
    label Nf = fuels.components().size();
    label Ns = sDB.composition().Y().size();

    // Calculate the interpolated gas properties at the position of the parcel
    vector Up = sDB.UInterpolator().interpolate(position(), cell())
        + Uturb();
    scalar rhog = sDB.rhoInterpolator().interpolate(position(), cell());
    scalar pg = sDB.pInterpolator().interpolate(position(), cell());
    scalar Tg = sDB.TInterpolator().interpolate(position(), cell());

    scalarField Yfg(Nf, 0.0);

    scalar cpMixture = 0.0;
    for(label i=0; i<Ns; i++)
    {
        const volScalarField& Yi = sDB.composition().Y()[i];
        if (sDB.isLiquidFuel()[i])
        {
            label j = sDB.gasToLiquidIndex()[i];
            scalar Yicelli = Yi[cell()];
            Yfg[j] = Yicelli;
        }
        cpMixture += Yi[cell()]*sDB.gasProperties()[i].Cp(Tg);
    }

    // correct the gaseous temperature for evaporated fuel

    scalar cellV = sDB.mesh().V()[cell()];
    scalar cellMass = rhog*cellV;
    Tg += sDB.shs()[cell()]/(cpMixture*cellMass);
    Tg = max(200.0, Tg);

    scalar tauMomentum = GREAT;
    scalar tauHeatTransfer = GREAT;
    scalarField tauEvaporation(Nf, GREAT);
    scalarField tauBoiling(Nf, GREAT);

    bool keepParcel = true;

    setRelaxationTimes
    (
        cell(),
        tauMomentum,
        tauEvaporation,
        tauHeatTransfer,
        tauBoiling,
        sDB,
        rhog,
        Up,
        Tg,
        pg,
        Yfg,
        m()*fuels.Y(X()),
        deltaT
    );


    // set the end-time for the track
    scalar tEnd = (1.0 - stepFraction())*deltaT;

    // set the maximum time step for this parcel
    scalar dtMax = min
    (
        tEnd,
        min
        (
            tauMomentum,
            min
            (
                1.0e+10*mag(min(tauEvaporation)), // evaporation is not an issue
                min
                (
                    mag(tauHeatTransfer),
                    mag(min(tauBoiling))
                )
            )
        )
    )/sDB.subCycles();

    // prevent the number of subcycles from being too many
    // (10 000 seems high enough)
    dtMax = max(dtMax, 1.0e-4*tEnd);

    bool switchProcessor = false;
    vector planeNormal = vector::zero;
    if (sDB.twoD())
    {
        planeNormal = n() ^ sDB.axisOfSymmetry();
        planeNormal /= mag(planeNormal);
    }

    // move the parcel until there is no 'timeLeft'
    while (keepParcel && tEnd > SMALL && !switchProcessor)
    {
        // set the lagrangian time-step
        scalar dt = min(dtMax, tEnd);

        // remember which cell the parcel is in
        // since this will change if a face is hit
        label celli = cell();
        scalar p = sDB.p()[celli];

        // track parcel to face, or end of trajectory
        if (keepParcel)
        {
            // Track and adjust the time step if the trajectory is not completed
            dt *= trackToFace(position() + dt*U_, sDB);

            // Decrement the end-time acording to how much time the track took
            tEnd -= dt;

            // Set the current time-step fraction.
            stepFraction() = 1.0 - tEnd/deltaT;

            if (onBoundary()) // hit face
            {
#               include "boundaryTreatment.H"
            }
        }

        if (keepParcel && sDB.twoD())
        {
            scalar z = position() & sDB.axisOfSymmetry();
            vector r = position() - z*sDB.axisOfSymmetry();
            if (mag(r) > SMALL)
            {
                correctNormal(sDB.axisOfSymmetry());
            }
        }

        // **** calculate the lagrangian source terms ****
        // First we get the 'old' properties.
        // and then 'update' them to get the 'new'
        // properties.
        // The difference is then added to the source terms.

        scalar oRho = fuels.rho(p, T(), X());
        scalarField oMass(Nf, 0.0);
        scalar oHg = 0.0;
        scalar oTotMass = m();
        scalarField oYf(fuels.Y(X()));

        forAll(oMass, i)
        {
            oMass[i] = m()*oYf[i];
            label j = sDB.liquidToGasIndex()[i];
            oHg += oYf[i]*sDB.gasProperties()[j].Hs(T());
        }

        vector oMom = m()*U();
        scalar oHv = fuels.hl(p, T(), X());
        scalar oH = oHg - oHv;
        scalar oPE = (p - fuels.pv(p, T(), X()))/oRho;

        // update the parcel properties (U, T, D)
        updateParcelProperties
        (
            dt,
            sDB,
            celli,
            face()
        );

        scalar nRho = fuels.rho(p, T(), X());
        scalar nHg = 0.0;
        scalarField nMass(Nf, 0.0);
        scalarField nYf(fuels.Y(X()));

        forAll(nMass, i)
        {
            nMass[i] = m()*nYf[i];
            label j = sDB.liquidToGasIndex()[i];
            nHg += nYf[i]*sDB.gasProperties()[j].Hs(T());
        }

        vector nMom = m()*U();
        scalar nHv = fuels.hl(p, T(), X());
        scalar nH = nHg - nHv;
        scalar nPE = (p - fuels.pv(p, T(), X()))/nRho;

        // Update the Spray Source Terms
        forAll(nMass, i)
        {
            sDB.srhos()[i][celli] += oMass[i] - nMass[i];
        }
        sDB.sms()[celli]   += oMom - nMom;

        sDB.shs()[celli]   +=
            oTotMass*(oH + oPE)
          - m()*(nH + nPE);

        // Remove evaporated mass from stripped mass
        ms() -= ms()*(oTotMass-m())/oTotMass;

        // remove parcel if it is 'small'
        if (m() < 1.0e-12)
        {
            keepParcel = false;

            // ... and add the removed 'stuff' to the gas
            forAll(nMass, i)
            {
                sDB.srhos()[i][celli] += nMass[i];
            }

            sDB.sms()[celli] += nMom;
            sDB.shs()[celli] += m()*(nH + nPE);
        }

        if (onBoundary() && keepParcel)
        {
            if (face() > -1)
            {
                if (isA<processorPolyPatch>(pbMesh[patch(face())]))
                {
                    switchProcessor = true;
                }
            }
        }
    }

    return keepParcel;
}


void Foam::parcel::updateParcelProperties
(
    const scalar dt,
    spray& sDB,
    const label celli,
    const label facei
)
{
    const liquidMixture& fuels = sDB.fuels();

    label Nf = sDB.fuels().components().size();
    label Ns = sDB.composition().Y().size();

    // calculate mean molecular weight
    scalar W = 0.0;
    for(label i=0; i<Ns; i++)
    {
        W += sDB.composition().Y()[i][celli]/sDB.gasProperties()[i].W();

    }
    W = 1.0/W;

    // Calculate the interpolated gas properties at the position of the parcel
    vector Up = sDB.UInterpolator().interpolate(position(), celli, facei)
        + Uturb();
    scalar rhog = sDB.rhoInterpolator().interpolate(position(), celli, facei);
    scalar pg = sDB.pInterpolator().interpolate(position(), celli, facei);
    scalar Tg0 = sDB.TInterpolator().interpolate(position(), celli, facei);

    // correct the gaseous temperature for evaporated fuel
    scalar cpMix = 0.0;
    for(label i=0; i<Ns; i++)
    {
        cpMix += sDB.composition().Y()[i][celli]
                *sDB.gasProperties()[i].Cp(T());
    }
    scalar cellV            = sDB.mesh().V()[celli];
    scalar rho              = sDB.rho()[celli];
    scalar cellMass         = rho*cellV;
    scalar dh               = sDB.shs()[celli];
    scalarField addedMass(Nf, 0.0);

    forAll(addedMass, i)
    {
        addedMass[i] += sDB.srhos()[i][celli]*cellV;
    }

    scalar Tg  = Tg0 + dh/(cpMix*cellMass);
    Tg = max(200.0, Tg);

    scalarField Yfg(Nf, 0.0);
    forAll(Yfg, i)
    {
        label j = sDB.liquidToGasIndex()[i];
        const volScalarField& Yj = sDB.composition().Y()[j];
        scalar Yfg0 = Yj[celli];
        Yfg[i] = (Yfg0*cellMass + addedMass[i])/(addedMass[i] + cellMass);
    }

    scalar tauMomentum     = GREAT;
    scalar tauHeatTransfer = GREAT;
    scalarField tauEvaporation(Nf, GREAT);
    scalarField tauBoiling(Nf, GREAT);

    setRelaxationTimes
    (
        celli,
        tauMomentum,
        tauEvaporation,
        tauHeatTransfer,
        tauBoiling,
        sDB,
        rhog,
        Up,
        Tg,
        pg,
        Yfg,
        m()*fuels.Y(X()),
        dt
    );

    scalar timeRatio = dt/tauMomentum;

    vector Ucorr = Up;
    vector gcorr = sDB.g();

    if (sDB.twoD())
    {
        // remove the tangential velocity component
        scalar v1 = Up & sDB.axisOfSymmetry();
        scalar v2 = Up & n();
        Ucorr     = v1*sDB.axisOfSymmetry() + v2*n();

        // Remove the tangential gravity component
        scalar g1 = gcorr & sDB.axisOfSymmetry();
        scalar g2 = gcorr & n();
        gcorr     = g1*sDB.axisOfSymmetry() + g2*n();
    }

    U() = (U() + timeRatio*Ucorr + gcorr*dt)/(1.0 + timeRatio);

    if (sDB.twoD())
    {
        vector normal = n() ^ sDB.axisOfSymmetry();
        normal /= mag(normal);
        scalar dU = U() & normal;
        U() -= dU*normal;
    }

    scalar TDroplet = T();
    scalar oldDensity = fuels.rho(pg, T(), X());
    scalar oldMass = m();
    scalarField Yf0(fuels.Y(X()));
    scalarField mi0(m()*Yf0);
    scalarField mi(mi0);

    scalar oldhg = 0.0;
    for (label i=0; i<Nf; i++)
    {
        label j = sDB.liquidToGasIndex()[i];
        oldhg += Yf0[i]*sDB.gasProperties()[j].Hs(T());
    }

    scalar oldhv = fuels.hl(pg, T(), X());
    scalar Np = N(oldDensity);

    scalar newDensity = oldDensity;
    scalar newMass    = oldMass;
    scalar newhg      = oldhg;
    scalar newhv      = oldhv;

    scalar Tnew = T();

    // NN.
    // first calculate the new temperature and droplet mass,
    // then calculate the energy source and correct the
    // gaseous temperature, Tg, and mass fraction, Yfg,
    // to calculate the new properties for the parcel
    // This procedure seems to be more stable
    label n = 0;
    while ((n < sDB.evaporation().nEvapIter()) && (m() > VSMALL))
    {
        n++;
        // new characteristic times does not need to be calculated the first time
        if (n > 1)
        {
            newDensity = fuels.rho(pg, Tnew, X());
            newMass = m();
            newhg = 0.0;
            scalarField Ynew(fuels.Y(X()));
            for(label i=0; i<Nf; i++)
            {
                label j = sDB.liquidToGasIndex()[i];
                newhg += Ynew[i]*sDB.gasProperties()[j].Hs(Tnew);
            }

            newhv = fuels.hl(pg, Tnew, X());

            scalar dm = oldMass - newMass;
            scalar dhNew = oldMass*(oldhg-oldhv) - newMass*(newhg-newhv);

            // Prediction of new gaseous temperature and fuel mass fraction
            Tg  = Tg0 + (dh+dhNew)/(cpMix*cellMass);
            Tg = max(200.0, Tg);

            forAll(Yfg, i)
            {
                label j = sDB.liquidToGasIndex()[i];
                const volScalarField& Yj = sDB.composition().Y()[j];
                scalar Yfg0 = Yj[celli];
                Yfg[i] = (Yfg0*cellMass + addedMass[i] + dm)
                        /(addedMass[i] + cellMass + dm);
            }

            setRelaxationTimes
            (
                celli,
                tauMomentum,
                tauEvaporation,
                tauHeatTransfer,
                tauBoiling,
                sDB,
                rhog,
                Up,
                Tg,
                pg,
                Yfg,
                m()*fuels.Y(X()),
                dt
            );

        }

        scalar Taverage = TDroplet + (Tg - TDroplet)/3.0;
        // for a liquid Cl \approx Cp
        scalar liquidcL = sDB.fuels().cp(pg, TDroplet, X());

        cpMix = 0.0;
        for(label i=0; i<Ns; i++)
        {
            if (sDB.isLiquidFuel()[i])
            {
                label j = sDB.gasToLiquidIndex()[i];
                cpMix += Yfg[j]*sDB.gasProperties()[i].Cp(Taverage);
            }
            else
            {
                scalar Y = sDB.composition().Y()[i][celli];
                cpMix += Y*sDB.gasProperties()[i].Cp(Taverage);
            }
        }

        scalar evaporationSource = 0.0;
        scalar z = 0.0;
        scalar tauEvap = 0.0;

        for(label i=0; i<Nf; i++)
        {
            tauEvap += X()[i]*fuels.properties()[i].W()/tauEvaporation[i];
        }
        tauEvap = fuels.W(X())/tauEvap;


        if (sDB.evaporation().evaporation())
        {
            scalar hv = fuels.hl(pg, TDroplet, X());
            evaporationSource =
                hv/liquidcL/tauEvap;

            z = cpMix*tauHeatTransfer/liquidcL/tauEvap;
        }

        if (sDB.heatTransfer().heatTransfer())
        {
            scalar fCorrect =
                sDB.heatTransfer().fCorrection(z)/tauHeatTransfer;

            Tnew =
                (TDroplet + dt*(fCorrect * Tg - evaporationSource))
                /(1.0 + dt*fCorrect);

            // Prevent droplet temperature to go above critial value
            Tnew = min(Tnew, fuels.Tc(X()));

            // Prevent droplet temperature to go too low
            // Mainly a numerical stability issue
            Tnew = max(200.0, Tnew);
            scalar Td = Tnew;

            scalar pAtSurface = fuels.pv(pg, Td, X());
            scalar pCompare = 0.999*pg;
            scalar boiling = pAtSurface >= pCompare;
            if (boiling)
            {
                // can not go above boiling temperature
                scalar Terr = 1.0e-3;
                label n=0;
                scalar dT = 1.0;
                scalar pOld = pAtSurface;
                while (dT > Terr)
                {
                    n++;
                    pAtSurface = fuels.pv(pg, Td, X());
                    if ((pAtSurface < pCompare) && (pOld < pCompare))
                    {
                        Td += dT;
                    }
                    else
                    {
                        if ((pAtSurface > pCompare) && (pOld > pCompare))
                        {
                            Td -= dT;
                        }
                        else
                        {
                            dT *= 0.5;
                            if ((pAtSurface > pCompare) && (pOld < pCompare))
                            {
                                Td -= dT;
                            }
                            else
                            {
                                Td += dT;
                            }
                        }
                    }
                    pOld = pAtSurface;
                    if (debug)
                    {
                        if (n>100)
                        {
                            Info << "n = " << n << ", T = " << Td << ", pv = " << pAtSurface << endl;
                        }
                    }
                }
                Tnew = Td;
            }
        }

        // Evaporate droplet!
        // if the droplet is NOT boiling use implicit scheme.
        if (sDB.evaporation().evaporation())
        {
            for(label i=0; i<Nf; i++)
            {
                // immediately evaporate mass that has reached critical
                // condition
                if (mag(Tnew - fuels.Tc(X())) < SMALL)
                {
                    mi[i] = 0.0;
                }
                else
                {
                    scalar Td = min(Tnew, 0.999*fuels.properties()[i].Tc());

                    scalar pAtSurface = fuels.properties()[i].pv(pg, Td);
                    scalar boiling = pAtSurface >= 0.999*pg;

                    if (!boiling)
                    {
                        scalar fr = dt/tauEvaporation[i];
                        mi[i] = mi0[i]/(1.0 + fr);
                    }
                    else
                    {
                        scalar fr = dt/tauBoiling[i];
                        mi[i] = mi0[i]/(1.0 + fr);
                    }
                }
            }

            scalar mTot = sum(mi);
            if (mTot > VSMALL)
            {
                scalarField Ynew(mi/mTot);
                scalarField Xnew(sDB.fuels().X(Ynew));
                forAll(Xnew, i)
                {
                    X()[i] = Xnew[i];
                }
                m() = mTot;
            }
            else
            {
                m() = 0.0;
            }
        }
        T() = Tnew;
        scalar rhod = fuels.rho(pg, T(), X());
        d() = cbrt(6.0*m()/(Np*rhod*M_PI));
    }

    T() = Tnew;
    scalar rhod = fuels.rho(pg, T(), X());
    m() = sum(mi);
    d() = cbrt(6.0*m()/(Np*rhod*M_PI));
}


void Foam::parcel::transformProperties(const tensor& T)
{
    U_ = transform(T, U_);
}


void Foam::parcel::transformProperties(const vector&)
{}


// ************************************************************************* //
