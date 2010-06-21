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

#include "hollowCone.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(hollowConeInjector, 0);

addToRunTimeSelectionTable
(
    injectorModel,
    hollowConeInjector,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
hollowConeInjector::hollowConeInjector
(
    const dictionary& dict,
    spray& sm
)
:
    injectorModel(dict, sm),
    hollowConeDict_(dict.subDict(typeName + "Coeffs")),
    dropletPDF_
    (
        pdfs::pdf::New
        (
            hollowConeDict_.subDict("dropletPDF"),
            sm.rndGen()
        )
    ),
    innerAngle_(hollowConeDict_.lookup("innerConeAngle")),
    outerAngle_(hollowConeDict_.lookup("outerConeAngle"))
{

    if (sm.injectors().size() != innerAngle_.size())
    {
        FatalError << "hollowConeInjector::hollowConeInjector"
            << "(const dictionary& dict, spray& sm)\n"
            << "Wrong number of entries in innerAngle"
            << abort(FatalError);
    }

    if (sm.injectors().size() != outerAngle_.size())
    {
        FatalError << "hollowConeInjector::hollowConeInjector"
            << "(const dictionary& dict, spray& sm)\n"
            << "Wrong number of entries in outerAngle"
            << abort(FatalError);
    }

    scalar referencePressure = sm.ambientPressure();

    // correct velocityProfile
    forAll(sm.injectors(), i)
    {
        sm.injectors()[i].properties()->correctProfiles(sm.fuels(), referencePressure);
    }

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

hollowConeInjector::~hollowConeInjector()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar hollowConeInjector::d0
(
    const label,
    const scalar
) const
{
    return dropletPDF_->sample();
}


vector hollowConeInjector::direction
(
    const label n,
    const label hole,
    const scalar time,
    const scalar d
) const
{
    scalar angle = innerAngle_[n] + rndGen_.scalar01()*(outerAngle_[n]-innerAngle_[n]);
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
    vector dir = dcorr*injectors_[n].properties()->direction(hole, time) + normal;
    dir /= mag(dir);

    return dir;

}

scalar hollowConeInjector::velocity
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

scalar hollowConeInjector::averageVelocity
(
    const label i
) const
{
    const injectorType& it = sm_.injectors()[i].properties();
    scalar dt = it.teoi() - it.tsoi();
    return it.integrateTable(it.velocityProfile())/dt;
}

} // End namespace Foam

// ************************************************************************* //
