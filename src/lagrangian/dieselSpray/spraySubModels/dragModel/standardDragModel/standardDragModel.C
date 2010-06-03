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

#include "standardDragModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(standardDragModel, 0);

addToRunTimeSelectionTable
(
    dragModel,
    standardDragModel,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
standardDragModel::standardDragModel
(
    const dictionary& dict
)
:
    dragModel(dict),
    dragDict_(dict.subDict(typeName + "Coeffs")),
    preReFactor_(readScalar(dragDict_.lookup("preReFactor"))),
    ReExponent_(readScalar(dragDict_.lookup("ReExponent"))),
    ReLimiter_(readScalar(dragDict_.lookup("ReLimiter"))),
    CdLimiter_(readScalar(dragDict_.lookup("CdLimiter"))),
    Cdistort_(readScalar(dragDict_.lookup("Cdistort")))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

standardDragModel::~standardDragModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar standardDragModel::Cd
(
    const scalar Re,
    const scalar dev
) const
{
    scalar drag = CdLimiter_;

    if (Re < ReLimiter_)
    {
        drag =  24.0*(1.0 + preReFactor_*pow(Re, ReExponent_))/Re;
    }

    // correct for deviation from sphericity
    drag *= (1.0 + Cdistort_*dev);

    return drag;

}


scalar standardDragModel::relaxationTime
(
    const vector& URel,
    const scalar diameter,
    const scalar rho,
    const scalar liquidDensity,
    const scalar nu,
    const scalar dev
) const
{

    scalar time = GREAT;
    scalar Re = mag(URel)*diameter/nu;

    if (Re > 0.1)
    {
        time = 4.0*liquidDensity*diameter /
        (
            3.0*rho*Cd(Re, dev)*mag(URel)
        );
    }
    else
    {
        // for small Re number, the relative velocity is both in
        // the nominator and denominator
        // use Cd = 24/Re and remove the SMALL/SMALL
        // expression for the velocities
        time = liquidDensity*diameter*diameter/(18*rho*nu*(1.0 + Cdistort_*dev));
    }
    return time;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
