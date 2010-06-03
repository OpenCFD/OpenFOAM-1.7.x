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

#include "multiHoleInjector.H"
#include "addToRunTimeSelectionTable.H"
#include "Random.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{

defineTypeNameAndDebug(multiHoleInjector, 0);

addToRunTimeSelectionTable
(
    injectorType,
    multiHoleInjector,
    dictionary
);
}
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::multiHoleInjector::multiHoleInjector
(
    const Foam::Time& t,
    const Foam::dictionary& dict
)
:
    injectorType(t, dict),
    propsDict_(dict.subDict(typeName + "Props")),
    centerPosition_(propsDict_.lookup("position")),
    xyAngle_(readScalar(propsDict_.lookup("xyAngle"))),
    zAngle_(readScalar(propsDict_.lookup("zAngle"))),
    nHoles_(readLabel(propsDict_.lookup("nHoles"))),
    umbrellaAngle_(readScalar(propsDict_.lookup("umbrellaAngle"))),
    nozzleTipDiameter_(readScalar(propsDict_.lookup("nozzleTipDiameter"))),
    angleSpacing_(propsDict_.lookup("angleSpacing")),
    d_(readScalar(propsDict_.lookup("diameter"))),
    Cd_(readScalar(propsDict_.lookup("Cd"))),
    mass_(readScalar(propsDict_.lookup("mass"))),
    nParcels_(readLabel(propsDict_.lookup("nParcels"))),
    X_(propsDict_.lookup("X")),
    massFlowRateProfile_(propsDict_.lookup("massFlowRateProfile")),
    velocityProfile_(massFlowRateProfile_),
    injectionPressureProfile_(massFlowRateProfile_),
    CdProfile_(massFlowRateProfile_),
    TProfile_(propsDict_.lookup("temperatureProfile")),
    averageParcelMass_(nHoles_*mass_/nParcels_),
    direction_(nHoles_),
    position_(nHoles_),
    pressureIndependentVelocity_(true),
    tangentialInjectionVector1_(nHoles_),
    tangentialInjectionVector2_(nHoles_)
{


    // check if time entries for soi and eoi match
    if (mag(massFlowRateProfile_[0][0]-TProfile_[0][0]) > SMALL)
    {
        FatalError << "multiHoleInjector::multiHoleInjector(const time& t, const dictionary dict) " << endl
            << " start-times do not match for TemperatureProfile and massFlowRateProfile."
            << abort(FatalError);
    }

    if (mag(massFlowRateProfile_[massFlowRateProfile_.size()-1][0]-TProfile_[TProfile_.size()-1][0]) > SMALL)
    {
        FatalError << "multiHoleInjector::multiHoleInjector(const time& t, const dictionary dict) " << endl
            << " end-times do not match for TemperatureProfile and massFlowRateProfile."
            << abort(FatalError);
    }

    // convert CA to real time
    forAll(massFlowRateProfile_, i)
    {
        massFlowRateProfile_[i][0] = t.userTimeToTime(massFlowRateProfile_[i][0]);
        velocityProfile_[i][0] = massFlowRateProfile_[i][0];
        injectionPressureProfile_[i][0] = massFlowRateProfile_[i][0];
    }
    
    forAll(TProfile_, i)
    {
        TProfile_[i][0] = t.userTimeToTime(TProfile_[i][0]);
    }

    scalar integratedMFR = integrateTable(massFlowRateProfile_);

    forAll(massFlowRateProfile_, i)
    {
        // correct the massFlowRateProfile to match the injected mass
        massFlowRateProfile_[i][1] *= mass_/integratedMFR;
        
        CdProfile_[i][0] = massFlowRateProfile_[i][0];
        CdProfile_[i][1] = Cd_;
    }

    setTangentialVectors();

    // check molar fractions
    scalar Xsum = 0.0;
    forAll(X_, i)
    {
        Xsum += X_[i];
    }

    if (mag(Xsum - 1.0) > SMALL)
    {
        Info << "Warning!!!\n multiHoleInjector::multiHoleInjector(const time& t, Istream& is)"
            << "X does not add up to 1.0, correcting molar fractions."
            << endl;
        forAll(X_, i)
        {
            X_[i] /= Xsum;
        }
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiHoleInjector::~multiHoleInjector()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiHoleInjector::setTangentialVectors()
{
    scalar pi180 = mathematicalConstant::pi/180.0;
    scalar alpha = xyAngle_*pi180;
    scalar phi = zAngle_*pi180;

    vector xp(cos(alpha), sin(alpha), 0.0);
    vector zp(cos(alpha)*sin(phi), sin(alpha)*sin(phi), cos(phi));
    if (mag(zp-xp) < 1.0e-15)
    {
        xp = vector(0.0, 0.0, -1.0);
        xp -= (xp & zp)*zp;
        xp /= mag(xp);
    }
    vector yp = zp^xp;

//    Info << "xp = " << xp << endl;
//    Info << "yp = " << yp << endl;
//    Info << "zp = " << zp << endl;

    scalar angle = 0.0;
    scalar u = umbrellaAngle_*pi180/2.0;
    for(label i=0; i<nHoles_; i++)
    {
        angle += angleSpacing_[i];
        scalar v = angle*pi180;
        direction_[i] = cos(v)*sin(u)*xp + sin(v)*sin(u)*yp + cos(u)*zp;
        vector dp = direction_[i] - (direction_[i] & zp)*direction_[i];
        if (mag(dp) > SMALL)
        {
            dp /= mag(dp);
        }
        position_[i] = centerPosition_ + 0.5*nozzleTipDiameter_*dp;
//        Info << "i = " << i << ", dir = " << direction_[i] << ", pos = " << position_[i] << endl;
    }

    Random rndGen(label(0));

    for(label i=0; i<nHoles_; i++)
    {
        vector tangent(vector::zero);
        scalar magV = 0;
        while (magV < SMALL)
        {
            vector testThis = rndGen.vector01();
            
            tangent = testThis - (testThis & direction_[i])*direction_[i];
            magV = mag(tangent);
        }

        tangentialInjectionVector1_[i] = tangent/magV;
        tangentialInjectionVector2_[i] = direction_[i] ^ tangentialInjectionVector1_[i];

    }   
}


Foam::label Foam::multiHoleInjector::nParcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{

    scalar mInj = mass_*(fractionOfInjection(time1)-fractionOfInjection(time0));
    label nParcels = label(mInj/averageParcelMass_ + 0.49);
    
    return nParcels;
}

const Foam::vector Foam::multiHoleInjector::position(const label n) const
{
    return position_[n];
}

Foam::vector Foam::multiHoleInjector::position
(
    const label n,
    const scalar time,
    const bool twoD,
    const scalar angleOfWedge,
    const vector& axisOfSymmetry,
    const vector& axisOfWedge,
    const vector& axisOfWedgeNormal,
    Random& rndGen
) const
{
    if (twoD)
    {
        scalar is = position_[n] & axisOfSymmetry;
        scalar magInj = mag(position_[n] - is*axisOfSymmetry);

        vector halfWedge =
            axisOfWedge*cos(0.5*angleOfWedge)
          + axisOfWedgeNormal*sin(0.5*angleOfWedge);
        halfWedge /= mag(halfWedge);

        return (is*axisOfSymmetry + magInj*halfWedge);
    }
    else
    {
        // otherwise, disc injection
        scalar iRadius = d_*rndGen.scalar01();
        scalar iAngle = 2.0*mathematicalConstant::pi*rndGen.scalar01();

        return
        ( 
            position_[n]
          + iRadius
          * (
              tangentialInjectionVector1_[n]*cos(iAngle)
            + tangentialInjectionVector2_[n]*sin(iAngle)
          )
        );
        
    }
    return position_[0];
}

Foam::label Foam::multiHoleInjector::nHoles() const
{
    return nHoles_;
}

Foam::scalar Foam::multiHoleInjector::d() const
{
    return d_;
}

const Foam::vector& Foam::multiHoleInjector::direction
(
    const label i,
    const scalar time
) const
{
    return direction_[i];
}

Foam::scalar Foam::multiHoleInjector::mass
(
    const scalar time0,
    const scalar time1,
    const bool twoD,
    const scalar angleOfWedge
) const
{
    scalar mInj = mass_*(fractionOfInjection(time1)-fractionOfInjection(time0));

    // correct mass if calculation is 2D 
    if (twoD)
    {
        mInj *= 0.5*angleOfWedge/mathematicalConstant::pi;
    }

    return mInj;
}

Foam::scalar Foam::multiHoleInjector::mass() const
{
    return mass_;
}

const Foam::scalarField& Foam::multiHoleInjector::X() const
{
    return X_;
}

Foam::List<Foam::multiHoleInjector::pair> Foam::multiHoleInjector::T() const
{
    return TProfile_;
}

Foam::scalar Foam::multiHoleInjector::T(const scalar time) const
{
    return getTableValue(TProfile_, time);
}

Foam::scalar Foam::multiHoleInjector::tsoi() const
{
    return massFlowRateProfile_[0][0];
}

Foam::scalar Foam::multiHoleInjector::teoi() const
{
    return massFlowRateProfile_[massFlowRateProfile_.size()-1][0];
}

Foam::scalar Foam::multiHoleInjector::massFlowRate
(
    const scalar time
) const
{
    return getTableValue(massFlowRateProfile_, time);
}

Foam::scalar Foam::multiHoleInjector::injectionPressure
(
    const scalar time
) const
{
    return getTableValue(injectionPressureProfile_, time);
}

Foam::scalar Foam::multiHoleInjector::velocity
(
    const scalar time
) const
{
    return getTableValue(velocityProfile_, time);
}

Foam::List<Foam::multiHoleInjector::pair> Foam::multiHoleInjector::CdProfile() const
{
    return CdProfile_;
}

Foam::scalar Foam::multiHoleInjector::Cd
(
    const scalar time
) const
{
    return Cd_;
}

Foam::scalar Foam::multiHoleInjector::fractionOfInjection(const scalar time) const
{
    return integrateTable(massFlowRateProfile_, time)/mass_;
}

Foam::scalar Foam::multiHoleInjector::injectedMass
(
    const scalar t
) const
{
    return mass_*fractionOfInjection(t);
}


void Foam::multiHoleInjector::correctProfiles
(
    const liquidMixture& fuel,
    const scalar referencePressure
)
{

    scalar A = nHoles_*0.25*mathematicalConstant::pi*pow(d_, 2.0);

    forAll(velocityProfile_, i)
    {
        scalar time = velocityProfile_[i][0];
        scalar rho = fuel.rho(referencePressure, T(time), X_);
        scalar v = massFlowRateProfile_[i][1]/(Cd_*rho*A);
        velocityProfile_[i][1] = v;
        injectionPressureProfile_[i][1] = referencePressure + 0.5*rho*v*v;
    }
}

Foam::vector Foam::multiHoleInjector::tan1(const label n) const
{
    return tangentialInjectionVector1_[n];
}

Foam::vector Foam::multiHoleInjector::tan2(const label n) const
{
    return tangentialInjectionVector2_[n];
}

// ************************************************************************* //
