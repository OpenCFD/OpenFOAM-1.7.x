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

#include "unitInjector.H"
#include "addToRunTimeSelectionTable.H"
#include "Random.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(unitInjector, 0);

    addToRunTimeSelectionTable
    (
        injectorType,
        unitInjector,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::unitInjector::unitInjector
(
    const Foam::Time& t,
    const Foam::dictionary& dict
)
:
    injectorType(t, dict),
    propsDict_(dict.subDict(typeName + "Props")),
    position_(propsDict_.lookup("position")),
    direction_(propsDict_.lookup("direction")),
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
    averageParcelMass_(mass_/nParcels_),
    pressureIndependentVelocity_(true)
{

    // check if time entries for soi and eoi match
    if (mag(massFlowRateProfile_[0][0]-TProfile_[0][0]) > SMALL)
    {
        FatalErrorIn
        (
            "unitInjector::unitInjector(const time& t, const dictionary dict)"
        )<< "start-times do not match for TemperatureProfile and "
         << " massFlowRateProfile." << nl << exit (FatalError);
    }

    if
    (
        mag(massFlowRateProfile_[massFlowRateProfile_.size()-1][0]
      - TProfile_[TProfile_.size()-1][0])
      > SMALL
    )
    {
        FatalErrorIn
        (
            "unitInjector::unitInjector(const time& t, const dictionary dict)"
        )<< "end-times do not match for TemperatureProfile and "
         << "massFlowRateProfile." << nl << exit(FatalError);
    }

    // convert CA to real time
    forAll(massFlowRateProfile_, i)
    {
        massFlowRateProfile_[i][0] =
            t.userTimeToTime(massFlowRateProfile_[i][0]);
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

    // Normalize the direction vector
    direction_ /= mag(direction_);

    setTangentialVectors();

    // check molar fractions
    scalar Xsum = 0.0;
    forAll(X_, i)
    {
        Xsum += X_[i];
    }

    if (mag(Xsum - 1.0) > SMALL)
    {
        WarningIn("unitInjector::unitInjector(const time& t, Istream& is)")
            << "X does not sum to 1.0, correcting molar fractions."
            << nl << endl;
        forAll(X_, i)
        {
            X_[i] /= Xsum;
        }
    }

}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::unitInjector::~unitInjector()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::unitInjector::setTangentialVectors()
{
    Random rndGen(label(0));
    scalar magV = 0.0;
    vector tangent;

    while (magV < SMALL)
    {
        vector testThis = rndGen.vector01();

        tangent = testThis - (testThis & direction_)*direction_;
        magV = mag(tangent);
    }

    tangentialInjectionVector1_ = tangent/magV;
    tangentialInjectionVector2_ = direction_ ^ tangentialInjectionVector1_;

}


Foam::label Foam::unitInjector::nParcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{
    scalar mInj = mass_*(fractionOfInjection(time1)-fractionOfInjection(time0));
    label nParcels = label(mInj/averageParcelMass_ + 0.49);
    return nParcels;
}


const Foam::vector Foam::unitInjector::position(const label n) const
{
    return position_;
}


Foam::vector Foam::unitInjector::position
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
        scalar is = position_ & axisOfSymmetry;
        scalar magInj = mag(position_ - is*axisOfSymmetry);

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
            position_
          + iRadius
          * (
              tangentialInjectionVector1_*cos(iAngle)
            + tangentialInjectionVector2_*sin(iAngle)
          )
        );

    }

    return position_;
}


Foam::label Foam::unitInjector::nHoles() const
{
    return 1;
}


Foam::scalar Foam::unitInjector::d() const
{
    return d_;
}


const Foam::vector& Foam::unitInjector::direction
(
    const label i,
    const scalar time
) const
{
    return direction_;
}


Foam::scalar Foam::unitInjector::mass
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


Foam::scalar Foam::unitInjector::mass() const
{
    return mass_;
}


const Foam::scalarField& Foam::unitInjector::X() const
{
    return X_;
}


Foam::List<Foam::unitInjector::pair> Foam::unitInjector::T() const
{
    return TProfile_;
}


Foam::scalar Foam::unitInjector::T(const scalar time) const
{
    return getTableValue(TProfile_, time);
}


Foam::scalar Foam::unitInjector::tsoi() const
{
    return massFlowRateProfile_[0][0];
}


Foam::scalar Foam::unitInjector::teoi() const
{
    return massFlowRateProfile_[massFlowRateProfile_.size()-1][0];
}


Foam::scalar Foam::unitInjector::massFlowRate(const scalar time) const
{
    return getTableValue(massFlowRateProfile_, time);
}


Foam::scalar Foam::unitInjector::injectionPressure(const scalar time) const
{
    return getTableValue(injectionPressureProfile_, time);
}


Foam::scalar Foam::unitInjector::velocity(const scalar time) const
{
    return getTableValue(velocityProfile_, time);
}


Foam::List<Foam::unitInjector::pair> Foam::unitInjector::CdProfile() const
{
    return CdProfile_;
}


Foam::scalar Foam::unitInjector::Cd(const scalar time) const
{
    return Cd_;
}


Foam::scalar Foam::unitInjector::fractionOfInjection(const scalar time) const
{
    return integrateTable(massFlowRateProfile_, time)/mass_;
}


Foam::scalar Foam::unitInjector::injectedMass(const scalar t) const
{
    return mass_*fractionOfInjection(t);
}


void Foam::unitInjector::correctProfiles
(
    const liquidMixture& fuel,
    const scalar referencePressure
)
{
    scalar A = 0.25*mathematicalConstant::pi*pow(d_, 2.0);
    scalar pDummy = 1.0e+5;

    forAll(velocityProfile_, i)
    {
        scalar time = velocityProfile_[i][0];
        scalar rho = fuel.rho(pDummy, T(time), X_);
        scalar v = massFlowRateProfile_[i][1]/(Cd_*rho*A);
        velocityProfile_[i][1] = v;
        injectionPressureProfile_[i][1] = referencePressure + 0.5*rho*v*v;
    }
}


Foam::vector Foam::unitInjector::tan1(const label) const
{
    return tangentialInjectionVector1_;
}


Foam::vector Foam::unitInjector::tan2(const label) const
{
    return tangentialInjectionVector2_;
}


// ************************************************************************* //
