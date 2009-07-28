/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "thermoPhysicsTypes.H"
#include "chemistrySolver.H"

#include "psiChemistryModel.H"
#include "rhoChemistryModel.H"

#include "EulerImplicit.H"
#include "ode.H"
#include "sequential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeChemistrySolver(psiChemistryModel, gasThermoPhysics)
    makeChemistrySolverType(EulerImplicit, psiChemistryModel, gasThermoPhysics)
    makeChemistrySolverType(ode, psiChemistryModel, gasThermoPhysics)
    makeChemistrySolverType(sequential, psiChemistryModel, gasThermoPhysics)

    makeChemistrySolver(psiChemistryModel, icoPoly8ThermoPhysics)
    makeChemistrySolverType
    (
        EulerImplicit,
        psiChemistryModel,
        icoPoly8ThermoPhysics
    )
    makeChemistrySolverType(ode, psiChemistryModel, icoPoly8ThermoPhysics)
    makeChemistrySolverType
    (
        sequential,
        psiChemistryModel,
        icoPoly8ThermoPhysics
    )

    makeChemistrySolver(rhoChemistryModel, gasThermoPhysics)
    makeChemistrySolverType(EulerImplicit, rhoChemistryModel, gasThermoPhysics)
    makeChemistrySolverType(ode, rhoChemistryModel, gasThermoPhysics)
    makeChemistrySolverType(sequential, rhoChemistryModel, gasThermoPhysics)

    makeChemistrySolver(rhoChemistryModel, icoPoly8ThermoPhysics)
    makeChemistrySolverType
    (
        EulerImplicit,
        rhoChemistryModel,
        icoPoly8ThermoPhysics
    )
    makeChemistrySolverType(ode, rhoChemistryModel, icoPoly8ThermoPhysics)
    makeChemistrySolverType
    (
        sequential,
        rhoChemistryModel,
        icoPoly8ThermoPhysics
    )
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
