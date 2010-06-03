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

#include "TAB.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(TAB, 0);

addToRunTimeSelectionTable
(
    breakupModel,
    TAB,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
TAB::TAB
(
    const dictionary& dict,
    spray& sm
)
:
    breakupModel(dict, sm),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    Cmu_(readScalar(coeffsDict_.lookup("Cmu"))),
    Comega_(readScalar(coeffsDict_.lookup("Comega"))),
    WeCrit_(readScalar(coeffsDict_.lookup("WeCrit")))
{

    // calculate the inverse function of the Rossin-Rammler Distribution
    const scalar xx0 = 12.0;
    const scalar rrd100 = 1.0/(1.0-exp(-xx0)*(1+xx0+pow(xx0, 2)/2+pow(xx0, 3)/6));

    for(label n=0; n<100; n++)
    {
        scalar xx = 0.12*(n+1);
        rrd_[n] = (1-exp(-xx)*(1 + xx + pow(xx, 2)/2 + pow(xx, 3)/6))*rrd100;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TAB::~TAB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void TAB::breakupParcel
(
    parcel& p,
    const scalar deltaT,
    const vector& Ug,
    const liquidMixture& fuels
) const
{

    scalar T = p.T();
    scalar pc = spray_.p()[p.cell()];
    scalar r = 0.5*p.d();
    scalar r2 = r*r;
    scalar r3 = r*r2;

    scalar rho = fuels.rho(pc, T, p.X());
    scalar sigma = fuels.sigma(pc, T, p.X());
    scalar mu = fuels.mu(pc, T, p.X());

    // inverse of characteristic viscous damping time
    scalar rtd = 0.5*Cmu_*mu/(rho*r2);
    
    // oscillation frequency (squared)
    scalar omega2 = Comega_*sigma/(rho*r3) - rtd*rtd;

    if (omega2 > 0)
    {
        scalar omega = sqrt(omega2);
        scalar rhog = spray_.rho()[p.cell()];
        scalar We = p.We(Ug, rhog, sigma);
        scalar Wetmp = We/WeCrit_;

        scalar y1 = p.dev() - Wetmp;
        scalar y2 = p.ddev()/omega;

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
                scalar coste = 1.0;
                if
                (
                    (Wetmp - a < -1)
                 && (p.ddev() < 0)
                )
                {
                    coste = -1.0;
                }
                
                scalar theta = acos((coste-Wetmp)/a);
                
                if (theta < phi)
                {
                    if (2*mathematicalConstant::pi-theta >= phi)
                    {
                        theta = -theta;
                    }
                    theta += 2*mathematicalConstant::pi;
                }
                tb = (theta-phi)/omega;

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
                scalar rs = r/
                (
                    1 
                  + (4.0/3.0)*pow(p.dev(), 2)
                  + rho*r3/(8*sigma)*pow(p.ddev(), 2)
                );
                
                label n = 0;
                bool found = false;
                scalar random = rndGen_.scalar01();
                while (!found && (n<99))
                {
                    if (rrd_[n]>random)
                    {
                        found = true;
                    }
                    n++;

                }
                scalar rNew = 0.04*n*rs;
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
