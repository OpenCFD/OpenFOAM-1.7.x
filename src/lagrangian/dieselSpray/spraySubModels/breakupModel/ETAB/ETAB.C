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

#include "ETAB.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ETAB, 0);

addToRunTimeSelectionTable
(
    breakupModel,
    ETAB,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ETAB::ETAB
(
    const dictionary& dict,
    spray& sm
)
:
    breakupModel(dict, sm),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    Cmu_(readScalar(coeffsDict_.lookup("Cmu"))),
    Comega_(readScalar(coeffsDict_.lookup("Comega"))),
    k1_(readScalar(coeffsDict_.lookup("k1"))),
    k2_(readScalar(coeffsDict_.lookup("k2"))),
    WeCrit_(readScalar(coeffsDict_.lookup("WeCrit"))),
    WeTransition_(readScalar(coeffsDict_.lookup("WeTransition"))),
    AWe_(0.0)
{
    scalar k21 = k2_/k1_;
    AWe_ = (k21*sqrt(WeTransition_) - 1.0)/pow(WeTransition_, 4.0);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ETAB::~ETAB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ETAB::breakupParcel
(
    parcel& p,
    const scalar deltaT,
    const vector& Ug,
    const liquidMixture& fuels
) const
{

    scalar T  = p.T();
    scalar pc  = spray_.p()[p.cell()];
    scalar r  = 0.5*p.d();
    scalar r2 = r*r;
    scalar r3 = r*r2;

    scalar rho   = fuels.rho(pc, T, p.X());
    scalar sigma = fuels.sigma(pc, T, p.X());
    scalar mu    = fuels.mu(pc, T, p.X());

    // inverse of characteristic viscous damping time
    scalar rtd = 0.5*Cmu_*mu/(rho*r2);

    // oscillation frequency (squared)
    scalar omega2 = Comega_*sigma/(rho*r3) - rtd*rtd;

    if (omega2 > 0)
    {
        scalar omega = sqrt(omega2);
        scalar romega = 1.0/omega;

        scalar rhog = spray_.rho()[p.cell()];
        scalar We = p.We(Ug, rhog, sigma);
        scalar Wetmp = We/WeCrit_;

        scalar y1 = p.dev() - Wetmp;
        scalar y2 = p.ddev()*romega;

        scalar a = sqrt(y1*y1 + y2*y2);

        // scotty we may have break-up
        if (a+Wetmp > 1.0)
        {
            scalar phic = y1/a;

            // constrain phic within -1 to 1
            phic = max(min(phic, 1), -1);

            scalar phit = acos(phic);
            scalar phi = phit;
            scalar quad = -y2/a;
            if (quad < 0)
            {
                phi = 2*mathematicalConstant::pi - phit;
            }

            scalar tb = 0;

            if (mag(p.dev()) < 1.0)
            {
                scalar theta = acos((1.0 - Wetmp)/a);

                if (theta < phi)
                {
                    if (2*mathematicalConstant::pi-theta >= phi)
                    {
                        theta = -theta;
                    }
                    theta += 2*mathematicalConstant::pi;
                }
                tb = (theta-phi)*romega;

                // breakup occurs
                if (deltaT > tb)
                {
                    p.dev() = 1.0;
                    p.ddev() = -a*omega*sin(omega*tb + phi);
                }
            }

            // update droplet size
            if (deltaT > tb)
            {
                scalar sqrtWe = AWe_*pow(We, 4.0) + 1.0;
                scalar Kbr = k1_*omega*sqrtWe;

                if (We > WeTransition_)
                {
                    sqrtWe = sqrt(We);
                    Kbr =k2_*omega*sqrtWe;
                }

                scalar rWetmp = 1.0/Wetmp;
                scalar cosdtbu = max(-1.0, min(1.0, 1.0-rWetmp));
                scalar dtbu = romega*acos(cosdtbu);
                scalar decay = exp(-Kbr*dtbu);

                scalar rNew = decay*r;
                if (rNew < r)
                {
                    p.d() = 2*rNew;
                    p.dev() = 0;
                    p.ddev() = 0;
                }
            }
        }
    }
    else
    {
        // reset droplet distortion parameters
        p.dev() = 0;
        p.ddev() = 0;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
