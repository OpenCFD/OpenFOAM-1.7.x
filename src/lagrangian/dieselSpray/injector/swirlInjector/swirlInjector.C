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

#include "swirlInjector.H"
#include "addToRunTimeSelectionTable.H"
#include "Random.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swirlInjector, 0);

    addToRunTimeSelectionTable
    (
        injectorType,
        swirlInjector,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swirlInjector::swirlInjector
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
    mass_(readScalar(propsDict_.lookup("mass"))),
    injectionPressure_(readScalar(propsDict_.lookup("injectionPressure"))),
    T_(readScalar(propsDict_.lookup("temperature"))),
    nParcels_(readLabel(propsDict_.lookup("nParcels"))),
    X_(propsDict_.lookup("X")),
    massFlowRateProfile_(propsDict_.lookup("massFlowRateProfile")),
    injectionPressureProfile_(propsDict_.lookup("injectionPressureProfile")),
    velocityProfile_(massFlowRateProfile_),
    CdProfile_(massFlowRateProfile_),
    TProfile_(massFlowRateProfile_),
    averageParcelMass_(mass_/nParcels_),
    pressureIndependentVelocity_(false)
{
    // convert CA to real time
    forAll(massFlowRateProfile_, i)
    {
        massFlowRateProfile_[i][0] =
            t.userTimeToTime(massFlowRateProfile_[i][0]);
    }
    forAll(injectionPressureProfile_, i)
    {
        injectionPressureProfile_[i][0] =
            t.userTimeToTime(injectionPressureProfile_[i][0]);
    }

    // check if time entries match
    if (mag(massFlowRateProfile_[0][0]-injectionPressureProfile_[0][0]) > SMALL)
    {
        FatalErrorIn
        (
            "swirlInjector::swirlInjector(const time& t, const dictionary dict)"
        )   << "Start-times do not match for "
               "injectionPressureProfile and massFlowRateProfile."
            << abort(FatalError);
    }

    // check if time entries match
    if
    (
        mag
        (
            massFlowRateProfile_[massFlowRateProfile_.size() - 1][0]
          - injectionPressureProfile_[injectionPressureProfile_.size() - 1][0]
        ) > SMALL
    )
    {
        FatalErrorIn
        (
            "swirlInjector::swirlInjector(const time& t, const dictionary dict)"
        )   << "End-times do not match for "
               "injectionPressureProfile and massFlowRateProfile."
            << abort(FatalError);
    }

    scalar integratedMFR = integrateTable(massFlowRateProfile_);
    scalar integratedPressure =
        integrateTable(injectionPressureProfile_)/(teoi()-tsoi());

    forAll(massFlowRateProfile_, i)
    {
        // correct the massFlowRateProfile to match the injected mass
        massFlowRateProfile_[i][1] *= mass_/integratedMFR;

        velocityProfile_[i][0] = massFlowRateProfile_[i][0];

        TProfile_[i][0] = massFlowRateProfile_[i][0];
        TProfile_[i][1] = T_;

        CdProfile_[i][0] = massFlowRateProfile_[i][0];

    }

    forAll(injectionPressureProfile_, i)
    {
        // correct the pressureProfile to match the injection pressure
        injectionPressureProfile_[i][1] *=
            injectionPressure_/integratedPressure;
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
        WarningIn
        (
            "swirlInjector::swirlInjector(const time& t, const dictionary dict)"
        )   << "X does not add up to 1.0, correcting molar fractions." << endl;

        forAll(X_, i)
        {
            X_[i] /= Xsum;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swirlInjector::~swirlInjector()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::swirlInjector::setTangentialVectors()
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

Foam::label Foam::swirlInjector::nParcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{

    scalar mInj = mass_*(fractionOfInjection(time1)-fractionOfInjection(time0));
    label nParcels = label(mInj/averageParcelMass_ + 0.49);

    return nParcels;
}

const Foam::vector Foam::swirlInjector::position(const label n) const
{
    return position_;
}

Foam::vector Foam::swirlInjector::position
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

Foam::label Foam::swirlInjector::nHoles() const
{
    return 1;
}

Foam::scalar Foam::swirlInjector::d() const
{
    return d_;
}

const Foam::vector& Foam::swirlInjector::direction
(
    const label i,
    const scalar time
) const
{
    return direction_;
}

Foam::scalar Foam::swirlInjector::mass
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

Foam::scalar Foam::swirlInjector::mass() const
{
    return mass_;
}

Foam::List<Foam::swirlInjector::pair>
Foam::swirlInjector::massFlowRateProfile() const
{
    return massFlowRateProfile_;
}

Foam::scalar Foam::swirlInjector::massFlowRate(const scalar time) const
{
    return getTableValue(massFlowRateProfile_, time);
}

Foam::List<Foam::swirlInjector::pair>
Foam::swirlInjector::injectionPressureProfile() const
{
    return injectionPressureProfile_;
}

Foam::scalar Foam::swirlInjector::injectionPressure(const scalar time) const
{
    return getTableValue(injectionPressureProfile_, time);
}

Foam::List<Foam::swirlInjector::pair>
Foam::swirlInjector::velocityProfile() const
{
    return velocityProfile_;
}

Foam::scalar Foam::swirlInjector::velocity(const scalar time) const
{
    return getTableValue(velocityProfile_, time);
}

Foam::List<Foam::swirlInjector::pair> Foam::swirlInjector::CdProfile() const
{
    return CdProfile_;
}

Foam::scalar Foam::swirlInjector::Cd(const scalar time) const
{
    return getTableValue(CdProfile_, time);
}

const Foam::scalarField& Foam::swirlInjector::X() const
{
    return X_;
}

Foam::List<Foam::swirlInjector::pair> Foam::swirlInjector::T() const
{
    return TProfile_;
}

Foam::scalar Foam::swirlInjector::T(const scalar time) const
{
    return T_;
}

Foam::scalar Foam::swirlInjector::tsoi() const
{
    return massFlowRateProfile_[0][0];
}

Foam::scalar Foam::swirlInjector::teoi() const
{
    return massFlowRateProfile_[massFlowRateProfile_.size()-1][0];
}

Foam::scalar Foam::swirlInjector::fractionOfInjection(const scalar time) const
{
    return integrateTable(massFlowRateProfile_, time)/mass_;
}

Foam::scalar Foam::swirlInjector::injectedMass
(
    const scalar t
) const
{
    return mass_*fractionOfInjection(t);
}

void Foam::swirlInjector::correctProfiles
(
    const liquidMixture& fuel,
    const scalar referencePressure
)
{

    scalar A = 0.25*mathematicalConstant::pi*pow(d_, 2.0);
    scalar pDummy = 1.0e+5;
    scalar rho = fuel.rho(pDummy, T_, X_);

    forAll(velocityProfile_, i)
    {
        scalar Pinj = getTableValue
        (
            injectionPressureProfile_,
            massFlowRateProfile_[i][0]
        );
        scalar mfr = massFlowRateProfile_[i][1]/(rho*A);
        scalar v = sqrt(2.0*(Pinj - referencePressure)/rho);
        velocityProfile_[i][1] = v;
        CdProfile_[i][1] = mfr/v;
    }
}

Foam::vector Foam::swirlInjector::tan1(const label n) const
{
    return tangentialInjectionVector1_;
}

Foam::vector Foam::swirlInjector::tan2(const label n) const
{
    return tangentialInjectionVector2_;
}


// ************************************************************************* //
