/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "parcel.H"

#include "spray.H"
#include "dragModel.H"
#include "evaporationModel.H"
#include "heatTransferModel.H"
#include "basicMultiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void parcel::setRelaxationTimes
(
    label celli,
    scalar& tauMomentum,
    scalarField& tauEvaporation,
    scalar& tauHeatTransfer,
    scalarField& tauBoiling,
    const spray& sDB,
    const scalar rho,
    const vector& Up,
    const scalar temperature,
    const scalar pressure,
    const scalarField& Yfg,
    const scalarField& m0,
    const scalar dt
)
{

    const liquidMixture& fuels = sDB.fuels();

    scalar mCell = rho*sDB.mesh().V()[cell()];
    scalarField mfg(Yfg*mCell);

    label Ns = sDB.composition().Y().size();
    label Nf = fuels.components().size();

    // Tf is based on the 1/3 rule
    scalar Tf  = T() + (temperature - T())/3.0;

    // calculate mixture properties
    scalar W = 0.0;
    scalar kMixture = 0.0;
    scalar cpMixture = 0.0;
    scalar muf = 0.0;

    for(label i=0; i<Ns; i++)
    {
        scalar Y = sDB.composition().Y()[i][celli];
        W += Y/sDB.gasProperties()[i].W();
        // Using mass-fractions to average...
        kMixture += Y*sDB.gasProperties()[i].kappa(Tf);
        cpMixture += Y*sDB.gasProperties()[i].Cp(Tf);
        muf += Y*sDB.gasProperties()[i].mu(Tf);
    }
    W = 1.0/W;

    scalarField Xf(Nf, 0.0);
    scalarField Yf(Nf, 0.0);
    scalarField psat(Nf, 0.0);
    scalarField msat(Nf, 0.0);

    for(label i=0; i<Nf; i++)
    {
        label j = sDB.liquidToGasIndex()[i];
        scalar Y = sDB.composition().Y()[j][celli];
        scalar Wi = sDB.gasProperties()[j].W();
        Yf[i] = Y;
        Xf[i] = Y*W/Wi;
        psat[i] = fuels.properties()[i].pv(pressure, temperature);
        msat[i] = min(1.0, psat[i]/pressure)*Wi/W;
    }

    scalar nuf = muf/rho;

    scalar liquidDensity = fuels.rho(pressure, T(), X());
    scalar liquidcL = fuels.cp(pressure, T(), X());
    scalar heatOfVapour = fuels.hl(pressure, T(), X());

    // calculate the partial rho of the fuel vapour
    // alternative is to use the mass fraction
    // however, if rhoFuelVap is small (zero)
    // d(mass)/dt = 0 => no evaporation... hmmm... is that good? NO!

    // Assume equilibrium at drop-surface => pressure @ surface
    // = vapour pressure to calculate fuel-vapour density @ surface
    scalar pressureAtSurface = fuels.pv(pressure, T(), X());
    scalar rhoFuelVap = pressureAtSurface*fuels.W(X())/(specie::RR*Tf);

    scalarField Xs(sDB.fuels().Xs(pressure, temperature, T(), Xf, X()));
    scalarField Ys(Nf, 0.0);
    scalar Wliq = 0.0;

    for(label i=0; i<Nf; i++)
    {
        label j = sDB.liquidToGasIndex()[i];
        scalar Wi = sDB.gasProperties()[j].W();
        Wliq += Xs[i]*Wi;
    }

    for(label i=0; i<Nf; i++)
    {
        label j = sDB.liquidToGasIndex()[i];
        scalar Wi = sDB.gasProperties()[j].W();
        Ys[i] = Xs[i]*Wi/Wliq;
    }

    scalar Reynolds = Re(Up, nuf);
    scalar Prandtl = Pr(cpMixture, muf, kMixture);

    // calculate the characteritic times

    if(liquidCore_> 0.5)
    {
//      no drag for parcels in the liquid core..
        tauMomentum = GREAT;
    }
    else
    {
        tauMomentum = sDB.drag().relaxationTime
        (
            Urel(Up),
            d(),
            rho,
            liquidDensity,
            nuf,
            dev()
        );
    }

    // store the relaxationTime since it is needed in some breakup models.
    tMom_ = tauMomentum;

    tauHeatTransfer = sDB.heatTransfer().relaxationTime
    (
        liquidDensity,
        d(),
        liquidcL,
        kMixture,
        Reynolds,
        Prandtl
    );

    // evaporation-properties are evaluated at averaged temperature
    // set the boiling conditions true if pressure @ surface is 99.9%
    // of the pressure
    // this is mainly to put a limit on the evaporation time,
    // since tauEvaporation is very very small close to the boiling point.

    for(label i=0; i<Nf; i++)
    {
        scalar Td = min(T(), 0.999*fuels.properties()[i].Tc());
        bool boiling = fuels.properties()[i].pv(pressure, Td) >= 0.999*pressure;
        scalar Di = fuels.properties()[i].D(pressure, Td);
        scalar Schmidt = Sc(nuf, Di);

        scalar partialPressure = Xf[i]*pressure;

//      saturated vapour
        if(partialPressure > psat[i])
        {
            tauEvaporation[i] = GREAT;
        }
//      not saturated vapour
        else
        {
            if (!boiling)
            {
                // for saturation evaporation, only use 99.99% for numerical robustness
                scalar dm = max(SMALL, 0.9999*msat[i] - mfg[i]);

                tauEvaporation[i] = sDB.evaporation().relaxationTime
                (
                    d(),
                    fuels.properties()[i].rho(pressure, Td),
                    rhoFuelVap,
                    Di,
                    Reynolds,
                    Schmidt,
                    Xs[i],
                    Xf[i],
                    m0[i],
                    dm,
                    dt
                );
            }
            else
            {
                scalar Nusselt =
                    sDB.heatTransfer().Nu(Reynolds, Prandtl);

//              calculating the boiling temperature of the liquid at ambient pressure
                scalar tBoilingSurface = Td;

                label Niter = 0;
                scalar deltaT = 10.0;
                scalar dp0 = fuels.properties()[i].pv(pressure, tBoilingSurface) - pressure;
                while ((Niter < 200) && (mag(deltaT) > 1.0e-3))
                {
                    Niter++;
                    scalar pBoil = fuels.properties()[i].pv(pressure, tBoilingSurface);
                    scalar dp = pBoil - pressure;
                    if ( (dp > 0.0) && (dp0 > 0.0) )
                    {
                        tBoilingSurface -= deltaT;
                    }
                    else
                    {
                        if ( (dp < 0.0) && (dp0 < 0.0) )
                        {
                            tBoilingSurface += deltaT;
                        }
                        else
                        {
                            deltaT *= 0.5;
                            if ( (dp > 0.0) && (dp0 < 0.0) )
                            {
                                tBoilingSurface -= deltaT;
                            }
                            else
                            {
                                tBoilingSurface += deltaT;
                            }
                        }
                    }
                    dp0 = dp;
                }

                scalar vapourSurfaceEnthalpy = 0.0;
                scalar vapourFarEnthalpy = 0.0;

                for(label k = 0; k < sDB.gasProperties().size(); k++)
                {
                    vapourSurfaceEnthalpy += sDB.composition().Y()[k][celli]*sDB.gasProperties()[k].H(tBoilingSurface);
                    vapourFarEnthalpy += sDB.composition().Y()[k][celli]*sDB.gasProperties()[k].H(temperature);
                }

                scalar kLiquid = fuels.properties()[i].K(pressure, 0.5*(tBoilingSurface+T()));

                tauBoiling[i] = sDB.evaporation().boilingTime
                (
                    fuels.properties()[i].rho(pressure, Td),
                    fuels.properties()[i].cp(pressure, Td),
                    heatOfVapour,
                    kMixture,
                    Nusselt,
                    temperature - T(),
                    d(),
                    liquidCore(),
                    sDB.runTime().value() - ct(),
                    Td,
                    tBoilingSurface,
                    vapourSurfaceEnthalpy,
                    vapourFarEnthalpy,
                    cpMixture,
                    temperature,
                    kLiquid
                );
            }

        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
