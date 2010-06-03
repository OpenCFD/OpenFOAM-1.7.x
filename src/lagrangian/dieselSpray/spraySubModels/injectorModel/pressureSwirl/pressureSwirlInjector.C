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

#include "pressureSwirlInjector.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pressureSwirlInjector, 0);

addToRunTimeSelectionTable
(
    injectorModel,
    pressureSwirlInjector,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
pressureSwirlInjector::pressureSwirlInjector
(
    const dictionary& dict,
    spray& sm
)
:
    injectorModel(dict, sm),
    pressureSwirlInjectorDict_(dict.subDict(typeName + "Coeffs")),

    coneAngle_(pressureSwirlInjectorDict_.lookup("ConeAngle")),
    coneInterval_(pressureSwirlInjectorDict_.lookup("ConeInterval")),
    maxKv_(pressureSwirlInjectorDict_.lookup("maxKv")),

    angle_(0.0)
{

    if (sm.injectors().size() != coneAngle_.size())
    {
        FatalError << "pressureSwirlInjector::pressureSwirlInjector"
            << "(const dictionary& dict, spray& sm)\n"
            << "Wrong number of entries in innerAngle"
            << abort(FatalError);
    }

    scalar referencePressure = sm.p().average().value();

    // correct velocityProfile
    forAll(sm.injectors(), i)
    {
        sm.injectors()[i].properties()->correctProfiles(sm.fuels(), referencePressure);
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

pressureSwirlInjector::~pressureSwirlInjector()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar pressureSwirlInjector::d0
(
    const label n, 
    const scalar t
) const
{
    const injectorType& it = injectors_[n].properties();

    scalar c = rndGen_.scalar01();
    angle_ = coneAngle_[n]  + 2.0 * coneInterval_[n] * (0.5 - c) ;

    angle_ *= mathematicalConstant::pi/360.0;

    scalar injectedMassFlow = it.massFlowRate(t);
    
    scalar cosAngle = cos(angle_);   

    scalar rhoFuel = sm_.fuels().rho(sm_.ambientPressure(), it.T(t), it.X()); 
    scalar injectorDiameter = it.d();  
     
    scalar deltaPressure = deltaPressureInj(t,n);
    scalar kV = kv(n, injectedMassFlow, deltaPressure);
    scalar v = kV * sqrt(2.0*deltaPressure/rhoFuel);    

    u_ = v * cosAngle;
    
    scalar A = injectedMassFlow/(mathematicalConstant::pi*rhoFuel*u_);

    return (injectorDiameter-sqrt(pow(injectorDiameter,2)-4.0*A))/2.0;
}

vector pressureSwirlInjector::direction
(
    const label n,
    const label hole,
    const scalar time,
    const scalar d
) const
{

    scalar alpha = sin(angle_);
    scalar dcorr = cos(angle_);
    scalar beta = 2.0*mathematicalConstant::pi*rndGen_.scalar01();

    // randomly distributed vector normal to the injection vector
    vector normal = vector::zero;
    
    if (sm_.twoD())
    {
        scalar reduce = 0.01;
        // correct beta if this is a 2D run
        // map it onto the 'angleOfWedge'

        beta *= (1.0-2.0*reduce)*sm_.angleOfWedge()/(2.0*mathematicalConstant::pi);
        beta += reduce*sm_.angleOfWedge();
        normal = alpha*
        (
            sm_.axisOfWedge()*cos(beta) +
            sm_.axisOfWedgeNormal()*sin(beta)
        );
    }
    else
    {
        normal = alpha*
        (
            injectors_[n].properties()->tan1(hole)*cos(beta) +
            injectors_[n].properties()->tan2(hole)*sin(beta)
        );
    }
    
    // set the direction of injection by adding the normal vector
    vector dir = dcorr*injectors_[n].properties()->direction(hole, time) + normal;
    dir /= mag(dir);

    return dir;
}


scalar pressureSwirlInjector::velocity
(
    const label i,
    const scalar time
) const
{
    return u_*sqrt(1.0 + pow(tan(angle_),2.0));
}

scalar pressureSwirlInjector::averageVelocity
(
    const label i
) const
{    

    const injectorType& it = sm_.injectors()[i].properties();

    scalar dt = it.teoi() - it.tsoi();

    scalar injectedMassFlow = it.mass()/(it.teoi()-it.tsoi());

    scalar injectionPressure = averagePressure(i);

    scalar Tav = it.integrateTable(it.T())/dt;
    scalar rhoFuel = sm_.fuels().rho(sm_.ambientPressure(), Tav, it.X());  

    scalar kV = kv(i, injectedMassFlow, injectionPressure);

    return  kV*sqrt(2.0*(injectionPressure-sm_.ambientPressure())/rhoFuel);
}


scalar pressureSwirlInjector::kv
(
    const label inj,
    const scalar massFlow,
    const scalar dPressure
) const
{

    const injectorType& it = injectors_[inj].properties();

    scalar coneAngle = coneAngle_[inj];

    coneAngle *= mathematicalConstant::pi/360.0;

    scalar cosAngle = cos(coneAngle);
    scalar Tav = it.integrateTable(it.T())/(it.teoi()-it.tsoi());

    scalar rhoFuel = sm_.fuels().rho(sm_.ambientPressure(), Tav, it.X()); 
    scalar injectorDiameter = it.d();  
     
    scalar kv = max
    (
        maxKv_[inj], 
        4.0*massFlow
        *
        sqrt(rhoFuel/2.0/dPressure)
        /
        (mathematicalConstant::pi*pow(injectorDiameter, 2.0)*rhoFuel*cosAngle)
    );

    return min(1.0, kv);   
}




scalar pressureSwirlInjector::deltaPressureInj(const scalar time, const label inj) const
{
    return injectors_[inj].properties()->injectionPressure(time) - sm_.ambientPressure();   
}

scalar pressureSwirlInjector::averagePressure(const label inj) const
{

    const injectorType& it = sm_.injectors()[inj].properties();

    scalar dt = it.teoi() - it.tsoi();
    return it.integrateTable(it.injectionPressureProfile())/dt;
}

} // End namespace Foam

// ************************************************************************* //
