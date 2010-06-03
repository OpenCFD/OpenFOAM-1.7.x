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

#include "error.H"
#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(breakupModel, 0);

defineRunTimeSelectionTable(breakupModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
breakupModel::breakupModel
(
    const dictionary& dict,
    spray& sm
)
:
    dict_(dict),
    spray_(sm),
    rndGen_(sm.rndGen()),
    includeOscillation_(dict_.lookup("includeOscillation")),
    TABcoeffsDict_(dict.subDict("TABCoeffs")),
    y0_(0.0),
    yDot0_(0.0),
    TABComega_(0.0),
    TABCmu_(0.0),
    TABWeCrit_(0.0)
{
    if (includeOscillation_)
    {
        y0_ = readScalar(TABcoeffsDict_.lookup("y0"));
        yDot0_ = readScalar(TABcoeffsDict_.lookup("yDot0"));
        TABComega_ = readScalar(TABcoeffsDict_.lookup("Comega"));
        TABCmu_ = readScalar(TABcoeffsDict_.lookup("Cmu"));
        TABWeCrit_ = readScalar(TABcoeffsDict_.lookup("WeCrit"));
    }
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

breakupModel::~breakupModel()
{}

// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void breakupModel::updateParcelProperties
(
    parcel& p,
    const scalar deltaT,
    const vector& Ug,
    const liquidMixture& fuels
) const
{

    if(includeOscillation_)
    {
    
        scalar T = p.T();
        scalar pc = spray_.p()[p.cell()];
        scalar r = 0.5 * p.d();
        scalar r2 = r*r;
        scalar r3 = r*r2;
    
        scalar rho = fuels.rho(pc, T, p.X());
        scalar sigma = fuels.sigma(pc, T, p.X());
        scalar mu = fuels.mu(pc, T, p.X());
    
        // inverse of characteristic viscous damping time    
        scalar rtd = 0.5*TABCmu_*mu/(rho*r2);
        
        // oscillation frequency (squared)
        scalar omega2 = TABComega_ * sigma /(rho*r3) - rtd*rtd;
        
        if(omega2 > 0)
        {

            scalar omega = sqrt(omega2);
            scalar rhog = spray_.rho()[p.cell()];
            scalar We = p.We(Ug, rhog, sigma);
            scalar Wetmp = We/TABWeCrit_;

            scalar y1 = p.dev() - Wetmp;
            scalar y2 = p.ddev()/omega;
                       
            // update distortion parameters
            scalar c = cos(omega*deltaT);
            scalar s = sin(omega*deltaT);
            scalar e = exp(-rtd*deltaT);
            y2 = (p.ddev() + y1*rtd)/omega;
            
            p.dev() = Wetmp + e*(y1*c + y2*s);
            if (p.dev() < 0)
            {
                p.dev() = 0.0;
                p.ddev() = 0.0;
            }
            else 
            {
                p.ddev() = (Wetmp-p.dev())*rtd + e*omega*(y2*c - y1*s);
            }
        }
        else
        {
            // reset droplet distortion parameters
            p.dev() = 0;
            p.ddev() = 0;
        }

    }
    
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
