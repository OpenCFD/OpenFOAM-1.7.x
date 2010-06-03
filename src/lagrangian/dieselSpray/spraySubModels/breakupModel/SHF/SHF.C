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

#include "SHF.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(SHF, 0);

addToRunTimeSelectionTable
(
    breakupModel,
    SHF,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
SHF::SHF
(
    const dictionary& dict,
    spray& sm
)
:
    breakupModel(dict, sm),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    g_(sm.g()),
    rndGen_(sm.rndGen()),
    weCorrCoeff_(readScalar(coeffsDict_.lookup("weCorrCoeff"))),
    weBuCrit_(readScalar(coeffsDict_.lookup("weBuCrit"))),
    weBuBag_(readScalar(coeffsDict_.lookup("weBuBag"))),
    weBuMM_(readScalar(coeffsDict_.lookup("weBuMM"))),
    ohnCoeffCrit_(readScalar(coeffsDict_.lookup("ohnCoeffCrit"))),
    ohnCoeffBag_(readScalar(coeffsDict_.lookup("ohnCoeffBag"))),
    ohnCoeffMM_(readScalar(coeffsDict_.lookup("ohnCoeffMM"))),
    ohnExpCrit_(readScalar(coeffsDict_.lookup("ohnExpCrit"))),
    ohnExpBag_(readScalar(coeffsDict_.lookup("ohnExpBag"))),
    ohnExpMM_(readScalar(coeffsDict_.lookup("ohnExpMM"))),
    cInit_(readScalar(coeffsDict_.lookup("Cinit"))),
    c1_(readScalar(coeffsDict_.lookup("C1"))),
    c2_(readScalar(coeffsDict_.lookup("C2"))),
    c3_(readScalar(coeffsDict_.lookup("C3"))),
    cExp1_(readScalar(coeffsDict_.lookup("Cexp1"))),
    cExp2_(readScalar(coeffsDict_.lookup("Cexp2"))),
    cExp3_(readScalar(coeffsDict_.lookup("Cexp3"))),
    weConst_(readScalar(coeffsDict_.lookup("Weconst"))),
    weCrit1_(readScalar(coeffsDict_.lookup("Wecrit1"))),
    weCrit2_(readScalar(coeffsDict_.lookup("Wecrit2"))),
    coeffD_(readScalar(coeffsDict_.lookup("CoeffD"))),
    onExpD_(readScalar(coeffsDict_.lookup("OnExpD"))),
    weExpD_(readScalar(coeffsDict_.lookup("WeExpD"))),
    mu_(readScalar(coeffsDict_.lookup("mu"))),
    sigma_(readScalar(coeffsDict_.lookup("sigma"))),
    d32Coeff_(readScalar(coeffsDict_.lookup("d32Coeff"))),
    cDmaxBM_(readScalar(coeffsDict_.lookup("cDmaxBM"))),
    cDmaxS_(readScalar(coeffsDict_.lookup("cDmaxS"))),
    corePerc_(readScalar(coeffsDict_.lookup("corePerc")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

SHF::~SHF()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void SHF::breakupParcel
(
    parcel& p,
    const scalar deltaT,
    const vector& vel,
    const liquidMixture& fuels
) const
{

    label celli = p.cell();
    scalar T = p.T();
    scalar pc = spray_.p()[celli];

    scalar sigma = fuels.sigma(pc, T, p.X());
    scalar rhoLiquid = fuels.rho(pc, T, p.X());
    scalar muLiquid = fuels.mu(pc, T, p.X());
    scalar rhoGas = spray_.rho()[celli];

    scalar weGas      = p.We(vel, rhoGas, sigma);
    scalar weLiquid   = p.We(vel, rhoLiquid, sigma);

    // correct the Reynolds number. Reitz is using radius instead of diameter

    scalar reLiquid   = p.Re(rhoLiquid, vel, muLiquid);
    scalar ohnesorge  = sqrt(weLiquid)/(reLiquid + VSMALL);

    vector acceleration = p.Urel(vel)/p.tMom();
    vector trajectory = p.U()/mag(p.U());

    vector vRel = p.Urel(vel);

    scalar weGasCorr = weGas/(1.0 + weCorrCoeff_ * ohnesorge);

    // droplet deformation characteristic time

    scalar tChar = p.d()/mag(vRel)*sqrt(rhoLiquid/rhoGas);

    scalar tFirst = cInit_ * tChar;

    scalar tSecond = 0;
    scalar tCharSecond = 0;

    bool bag = false;
    bool multimode = false;
    bool shear = false;
    bool success = false;


    //  updating the droplet characteristic time
    p.ct() += deltaT;

    if(weGas > weConst_)
    {
        if(weGas < weCrit1_)
        {
            tCharSecond = c1_*pow((weGas - weConst_),cExp1_);
        }
        else if(weGas >= weCrit1_ && weGas <= weCrit2_)
        {
            tCharSecond = c2_*pow((weGas - weConst_),cExp2_);
        }
        else
        {
            tCharSecond = c3_*pow((weGas - weConst_),cExp3_);
        }
    }

    scalar weC = weBuCrit_*(1.0+ohnCoeffCrit_*pow(ohnesorge,ohnExpCrit_));
    scalar weB = weBuBag_*(1.0+ohnCoeffBag_*pow(ohnesorge, ohnExpBag_));
    scalar weMM = weBuMM_*(1.0+ohnCoeffMM_*pow(ohnesorge, ohnExpMM_));

    if(weGas > weC && weGas < weB)
    {
        bag = true;
    }

    if(weGas >= weB && weGas <= weMM)
    {
        multimode = true;
    }

    if(weGas > weMM)
    {
        shear = true;
    }

    tSecond = tCharSecond * tChar;

    scalar tBreakUP = tFirst + tSecond;
    if(p.ct() > tBreakUP)
    {

        scalar d32 = coeffD_*p.d()*pow(ohnesorge,onExpD_)*pow(weGasCorr,weExpD_);

        if(bag || multimode)
        {

            scalar d05 = d32Coeff_ * d32;

            scalar x = 0.0;
            scalar y = 0.0;
            scalar d = 0.0;

            while(!success)
            {
                x = cDmaxBM_*rndGen_.scalar01();
                d = sqr(x)*d05;
                y = rndGen_.scalar01();

                scalar p = x/(2.0*sqrt(2.0*mathematicalConstant::pi)*sigma_)*exp(-0.5*sqr((x-mu_)/sigma_));

                if (y<p)
                {
                    success = true;
                }
            }

            p.d() = d;
            p.ct() = 0.0;
        }

        if(shear)
        {
            scalar dC = weConst_*sigma/(rhoGas*sqr(mag(vRel)));
            scalar d32Red = 4.0*(d32 * dC)/(5.0 * dC - d32);
            scalar initMass = p.m();

            scalar d05 = d32Coeff_ * d32Red;

            scalar x = 0.0;
            scalar y = 0.0;
            scalar d = 0.0;

            while(!success)
            {

                x = cDmaxS_*rndGen_.scalar01();
                d = sqr(x)*d05;
                y = rndGen_.scalar01();

                scalar p = x/(2.0*sqrt(2.0*mathematicalConstant::pi)*sigma_)*exp(-0.5*sqr((x-mu_)/sigma_));

                if (y<p)
                {
                    success = true;
                }
            }

            p.d() = dC;
            p.m() = corePerc_ * initMass;

            spray_.addParticle
            (
                new parcel
                (
                    spray_,
                    p.position(),
                    p.cell(),
                    p.n(),
                    d,
                    p.T(),
                    (1.0 - corePerc_)* initMass,
                    0.0,
                    0.0,
                    0.0,
                    -GREAT,
                    p.tTurb(),
                    0.0,
                    scalar(p.injector()),
                    p.U(),
                    p.Uturb(),
                    p.X(),
                    p.fuelNames()
                )
            );

            p.ct() = 0.0;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
