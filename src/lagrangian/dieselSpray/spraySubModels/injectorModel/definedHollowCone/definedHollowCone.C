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

#include "definedHollowCone.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(definedHollowConeInjector, 0);

addToRunTimeSelectionTable
(
    injectorModel,
    definedHollowConeInjector,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
definedHollowConeInjector::definedHollowConeInjector
(
    const dictionary& dict,
    spray& sm
)
:
    injectorModel(dict, sm),
    definedHollowConeDict_(dict.subDict(typeName + "Coeffs")),
    dropletPDF_
    (
        pdfs::pdf::New
        (
            definedHollowConeDict_.subDict("dropletPDF"),
            sm.rndGen()
        )
    ),
    innerConeAngle_(definedHollowConeDict_.lookup("innerConeAngle")),
    outerConeAngle_(definedHollowConeDict_.lookup("outerConeAngle"))
{

    // convert CA to real time - inner cone angle
    forAll(innerConeAngle_, i)
    {
        innerConeAngle_[i][0] = sm.runTime().userTimeToTime(innerConeAngle_[i][0]);
    }
    // convert CA to real time - outer cone angle
    forAll(outerConeAngle_, i)
    {
        outerConeAngle_[i][0] = sm.runTime().userTimeToTime(outerConeAngle_[i][0]);
    }

    // check number of injectors
    if (sm.injectors().size() != 1)
    {
        Info << "Warning!!!\n"
             << "definedHollowConeInjector::definedHollowConeInjector"
             << "(const dictionary& dict, spray& sm)\n"
             << "Same inner/outer cone angle profiles applied to each injector"
             << endl;
    }

    // check number of entries in innerConeAngle list
    if (innerConeAngle_.empty())
    {
        FatalError << "definedHollowConeInjector::definedHollowConeInjector"
             << "(const dictionary& dict, spray& sm)\n"
             << "Number of entries in innerConeAngle must be greater than zero"
             << abort(FatalError);
    }

    // check number of entries in outerConeAngle list
    if (outerConeAngle_.empty())
    {
        FatalError << "definedHollowConeInjector::definedHollowConeInjector"
             << "(const dictionary& dict, spray& sm)\n"
             << "Number of entries in outerConeAngle must be greater than zero"
             << abort(FatalError);
    }

    scalar referencePressure = sm.p().average().value();
    // correct pressureProfile
    forAll(sm.injectors(), i)
    {
        sm.injectors()[i].properties()->correctProfiles(sm.fuels(), referencePressure);
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

definedHollowConeInjector::~definedHollowConeInjector()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar definedHollowConeInjector::d0
(
    const label n,
    const scalar t
) const
{
    // swallow function arguments - not used
    // return value sampled from PDF
    return dropletPDF_->sample();
}


vector definedHollowConeInjector::direction
(
    const label n,
    const label hole,
    const scalar t,
    const scalar d
) const
{

    const injectorType& it = injectors_[n].properties();

    // interpolate to find inner and outer angles at time, t
    scalar angleInner = it.getTableValue(innerConeAngle_, t);
    scalar angleOuter = it.getTableValue(outerConeAngle_, t);

    // use random number to generate angle between inner/outer cone angles
    scalar angle = angleInner + rndGen_.scalar01()*(angleOuter-angleInner);

    scalar alpha = sin(angle*mathematicalConstant::pi/360.0);
    scalar dcorr = cos(angle*mathematicalConstant::pi/360.0);
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
    vector dir = dcorr*injectors_[n].properties()->direction(hole, t) + normal;
    // normailse direction vector
    dir /= mag(dir);

    return dir;

}

scalar definedHollowConeInjector::velocity
(
    const label i,
    const scalar time
) const
{
    const injectorType& it = sm_.injectors()[i].properties();
    if (it.pressureIndependentVelocity())
    {
        return it.getTableValue(it.velocityProfile(), time);
    }
    else
    {
        scalar Pref = sm_.ambientPressure();
        scalar Pinj = it.getTableValue(it.injectionPressureProfile(), time);
        scalar rho = sm_.fuels().rho(Pinj, it.T(time), it.X());
        scalar dp = max(0.0, Pinj - Pref);
        return sqrt(2.0*dp/rho);
    }
}

scalar definedHollowConeInjector::averageVelocity
(
    const label i
) const
{
    const injectorType& it = sm_.injectors()[i].properties();
    scalar dt = it.teoi() - it.tsoi();
    return it.integrateTable(it.velocityProfile())/dt;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
