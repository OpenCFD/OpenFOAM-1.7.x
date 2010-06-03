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

#include "saturateEvaporationModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(saturateEvaporationModel, 0);

addToRunTimeSelectionTable
(
    evaporationModel,
    saturateEvaporationModel,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
saturateEvaporationModel::saturateEvaporationModel
(
    const dictionary& dict
)
:
    evaporationModel(dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

saturateEvaporationModel::~saturateEvaporationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool saturateEvaporationModel::evaporation() const
{
    return true;
}

// Correlation for the Sherwood Number
scalar saturateEvaporationModel::Sh
(
    const scalar ReynoldsNumber,
    const scalar SchmidtNumber
) const
{
    return 0.0;
}

scalar saturateEvaporationModel::relaxationTime
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
    return max(SMALL,dt*(m0/dm - 1.0));
}


scalar saturateEvaporationModel::boilingTime
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

inline label saturateEvaporationModel::nEvapIter() const
{
    return 1;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
