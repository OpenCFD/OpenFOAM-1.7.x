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

Class
    Foam::BasicReactingCloud

Description
    Reacting cloud templated on the type of carrier phase thermodynamics

SourceFiles
    BasicReactingCloud.C

\*---------------------------------------------------------------------------*/

#ifndef BasicReactingCloud_H
#define BasicReactingCloud_H

#include "ReactingCloud.H"
#include "BasicReactingParcel.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef ReactingCloud<BasicReactingParcel<constGasThermoPhysics> >
        constThermoReactingCloud;

    typedef ReactingCloud<BasicReactingParcel<gasThermoPhysics> >
        thermoReactingCloud;

    typedef ReactingCloud<BasicReactingParcel<icoPoly8ThermoPhysics> >
        icoPoly8ThermoReactingCloud;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
