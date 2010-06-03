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

#include "standardEvaporationModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardEvaporationModel, 0);

addToRunTimeSelectionTable
(
    evaporationModel,
    standardEvaporationModel,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
standardEvaporationModel::standardEvaporationModel
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

standardEvaporationModel::~standardEvaporationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool standardEvaporationModel::evaporation() const
{
    return true;
}

// Correlation for the Sherwood Number
scalar standardEvaporationModel::Sh
(
    const scalar ReynoldsNumber,
    const scalar SchmidtNumber
) const
{
    return 2.0 + preReScFactor_*pow(ReynoldsNumber,ReExponent_)*pow(SchmidtNumber,ScExponent_);
}

scalar standardEvaporationModel::relaxationTime
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

    scalar denominator =
        6.0 * massDiffusionCoefficient
      * Sh(ReynoldsNumber, SchmidtNumber)
      * rhoFuelVapor * lgExpr;

    if (denominator > SMALL)
    {
        time = max(VSMALL, liquidDensity * pow(diameter, 2.0)/denominator);
    }

    return time;
}


scalar standardEvaporationModel::boilingTime
(
    const scalar liquidDensity,
    const scalar cpFuel,
    const scalar heatOfVapour,
    const scalar kappa,
    const scalar Nusselt,
    const scalar deltaTemp,
    const scalar diameter,
    const scalar, 
    const scalar, 
    const scalar, 
    const scalar, 
    const scalar, 
    const scalar, 
    const scalar, 
    const scalar, 
    const scalar 
) const
{
    scalar time = GREAT;

    // the deltaTemperature is limited to not go below .5 deg K
    // for numerical reasons.
    // This is probably not important, since it results in an upper
    // limit for the boiling time... which we have anyway.
    scalar deltaT = max(0.5, deltaTemp);

    time = liquidDensity*cpFuel*sqr(diameter)/
    (
        6.0 * kappa * Nusselt * log(1.0 + cpFuel * deltaT/max(SMALL, heatOfVapour))
    );

    time = max(VSMALL, time);

    return time;
}

inline label standardEvaporationModel::nEvapIter() const
{
    return nEvapIter_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
