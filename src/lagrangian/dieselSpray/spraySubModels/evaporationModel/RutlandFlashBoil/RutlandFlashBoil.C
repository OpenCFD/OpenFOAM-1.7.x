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

#include "RutlandFlashBoil.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(RutlandFlashBoil, 0);

addToRunTimeSelectionTable
(
    evaporationModel,
    RutlandFlashBoil,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
RutlandFlashBoil::RutlandFlashBoil
(
    const dictionary& dict
)
:
    evaporationModel(dict),
    evapDict_(dict.subDict(typeName + "Coeffs")),
    preReScFactor_(readScalar(evapDict_.lookup("preReScFactor"))),
    ReExponent_(readScalar(evapDict_.lookup("ReExponent"))),
    ScExponent_(readScalar(evapDict_.lookup("ScExponent"))),
    evaporationScheme_(evapDict_.lookup("evaporationScheme")),
    nEvapIter_(0)
{
    if (evaporationScheme_ == "implicit") 
    {
        nEvapIter_ = 2;
    }
    else if (evaporationScheme_ == "explicit") 
    {
        nEvapIter_ = 1;
    }
    else 
    {
        FatalError
            << "evaporationScheme type " << evaporationScheme_
            << " unknown.\n"
            << "Use implicit or explicit."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

RutlandFlashBoil::~RutlandFlashBoil()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool RutlandFlashBoil::evaporation() const
{
    return true;
}

// Correlation for the Sherwood Number
scalar RutlandFlashBoil::Sh
(
    const scalar ReynoldsNumber,
    const scalar SchmidtNumber
) const
{
    return 2.0 + preReScFactor_*pow(ReynoldsNumber,ReExponent_)*pow(SchmidtNumber,ScExponent_);
}

scalar RutlandFlashBoil::relaxationTime
(
    const scalar diameter,
    const scalar liquidDensity,
    const scalar rhoFuelVapor,
    const scalar massDiffusionCoefficient,
    const scalar ReynoldsNumber,
    const scalar SchmidtNumber,
    const scalar Xs,
    const scalar Xf,
    const scalar m0,
    const scalar dm,
    const scalar dt
) const
{
    scalar time = GREAT;
    scalar lgExpr = 0.0;

    /*
        (pressure - partialFuelVaporPressure)/
        (pressure - pressureAtSurface)
      = 1 + Xratio

        if the pressure @ Surface > pressure
        this lead to boiling
        and Xratio -> infinity (as it should)
        ... this is numerically nasty

    NB! by N. Nordin
        X_v,s = (p_v,s/p) X_v,d
        where X_v,d = 1 for single component fuel
        according to eq (3.136)
        in D. Clerides Thesis
    */

    scalar Xratio = (Xs - Xf)/max(SMALL, 1.0 - Xs);

    if (Xratio > 0.0)
    {
        lgExpr = log(1.0 + Xratio);
    }

    // From Equation (3.79) in C. Kralj's Thesis:
    // Note that the 2.0 (instead of 6.0) below is correct, since evaporation
    // is d(diameter)/dt and not d(mass)/dt
    
    scalar Sherwood = Sh(ReynoldsNumber, SchmidtNumber);
    
    scalar FbExp = 0.7;
    
    scalar logXratio = log(1.0+Xratio);
    scalar Fb = 1.0;
    
    if(logXratio > SMALL)
    {
        Fb = pow((1.0 + Xratio),FbExp) * log(1.0+Xratio)/Xratio;
    }

//  TL: proposed correction to sherwood number, implemented

    Sherwood = 2.0 + (Sherwood - 2.0)/Fb;

    scalar denominator =
        6.0 * massDiffusionCoefficient
      * Sherwood
      * rhoFuelVapor * lgExpr;

    if (denominator > SMALL)
    {
        time = max(VSMALL, liquidDensity * pow(diameter, 2.0)/denominator);
    }

    return time;
}


scalar RutlandFlashBoil::boilingTime
(
    const scalar liquidDensity,
    const scalar cpFuel,
    const scalar heatOfVapour,
    const scalar kappa,
    const scalar Nusselt,
    const scalar deltaTemp,
    const scalar diameter,
    const scalar liquidCore,
    const scalar ct,
    const scalar tDrop,
    const scalar tBoilingSurface,
    const scalar vapourSurfaceEnthalpy,
    const scalar vapourFarEnthalpy,
    const scalar cpGas,
    const scalar temperature,
    const scalar kLiq
) const
{

    scalar time = GREAT;

    // the deltaTemperature is limited to not go below .5 deg K
    // for numerical reasons.
    // This is probably not important, since it results in an upper
    // limit for the boiling time... which we have anyway.

    //  TL kSet to the k value at the droplet temperature, not as in the Rutland Paper
    
    if(liquidCore > 0.5)
    {
        if(tDrop > tBoilingSurface)
        {               
            //  Evaporation of the liquid sheet      
           
            scalar psi = 2.72;
            scalar kIncreased = psi * kLiq;
            scalar alfa = psi * kIncreased/(liquidDensity * cpFuel);
            scalar F = alfa * ct/sqr(0.5 * diameter);
    
            scalar expSum = 0.0;
            scalar expSumOld = expSum;
        
            label Niter = 200;
                
            for(label k=0; k < Niter ; k++)
            {
                expSum += exp(sqr(-k*mathematicalConstant::pi*sqrt(F)/2.0));
                if(mag(expSum-expSumOld)/expSum < 1.0e-3)
                {
                    break;    
                }
                expSumOld = expSum;
            }
        }
    }
    else
    {
        scalar dTLB = min(0.5, tDrop -  tBoilingSurface);
        scalar alfaS = 0.0;

        if(dTLB >= 0.0 &&  dTLB < 5.0)
        {
            alfaS = 0.76 * pow(dTLB, 0.26);
        } 
        if(dTLB >= 5.0 &&  dTLB < 25.0)
        {
            alfaS = 0.027 * pow(dTLB, 2.33);
        }
        if(dTLB >= 25.0)
        {
            alfaS = 13.8 * pow(dTLB, 0.39);
        } 
        
        scalar Gf = 
        (
            4.0 * alfaS * dTLB * mathematicalConstant::pi * sqr(diameter/2.0)
        ) 
        / 
        heatOfVapour; 
        
        //  calculation of the heat transfer vapourization at superheated conditions (temperature>tBoilingSurface)
        scalar G = 0.0;
        if(temperature > tBoilingSurface)
        {
            scalar NusseltCorr = Nusselt ;
            scalar A = mag((vapourFarEnthalpy-vapourSurfaceEnthalpy)/heatOfVapour);

            // TL : 2.0? or 1.0? try 1!
            scalar B = 1.0*mathematicalConstant::pi*kappa/cpGas*diameter*NusseltCorr;
            scalar nPos = B * log(1.0 + A)/Gf + 1.0;  
            scalar nNeg = (1.0/A)*(exp(Gf/B) - 1.0 - A) + 1.0;
        
            scalar Gpos = Gf*nPos;
            scalar Gneg = Gf/nNeg;
            
            //scalar FgPos = Gpos + Gf - B * log( 1.0 + ( 1.0 + Gf/Gpos ) * A);
            scalar FgNeg = Gneg + Gf - B * log( 1.0 + ( 1.0 + Gf/Gneg ) * A);
            
            if(FgNeg > 0.0)
            {
                for(label j = 0; j < 20; j++)
                {
                    Gneg = Gneg/10.0;
                    Gneg = max(Gneg, VSMALL);
                    FgNeg = Gneg + Gf - B * log( 1.0 + ( 1.0 + Gf/Gneg ) * A);
                    if(FgNeg < 0.0)
                    {
                        break;
                    }            
                }
            }        
            
            FgNeg = Gneg + Gf - B * log( 1.0 + ( 1.0 + Gf/Gneg ) * A);                
            
            G = 0.5*(Gpos+Gneg);
            scalar Gold = -100;
        
            label Niter = 200;
            label k=0;
            
            if(FgNeg > 0.0)
            {
                Info << "no convergence" << endl;   
            }
            
            
            if(FgNeg < 0.0)
            {
            
                for(k=0; k<Niter ; k++)
                {

                    scalar Fg = G + Gf - B * log( 1.0 + ( 1.0 + Gf/G ) * A);

                    if(Fg > 0)
                    {
                        Gpos = G;
                        G = 0.5*(Gpos+Gneg);
                    }
                    else
                    {
                        Gneg = G;
                        G = 0.5*(Gpos+Gneg);
                    } 
            
                    Gold = G;
                    if(mag(G-Gold)/Gold < 1.0e-3)
                    {
                        break;
                    }
                }
        
                if(k >= Niter - 1)
                {
                    Info << " No convergence for G " << endl;
                }
            }
            else
            {
                G = 0.0;
            }
        }
        
        time = ((4.0/3.0)*mathematicalConstant::pi*pow(diameter/2.0,3.0))*liquidDensity/(G+Gf);
        time = max(VSMALL, time);
    }

    return time;
}

inline label RutlandFlashBoil::nEvapIter() const
{
    return nEvapIter_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
