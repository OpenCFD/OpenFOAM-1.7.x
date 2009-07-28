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

#include "dimensionSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::dimensionSet Foam::dimless(0, 0, 0, 0, 0, 0, 0);

const Foam::dimensionSet Foam::dimMass(1, 0, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimLength(0, 1, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimTime(0, 0, 1, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimTemperature(0, 0, 0, 1, 0, 0, 0);
const Foam::dimensionSet Foam::dimMoles(0, 0, 0, 0, 1, 0, 0);
const Foam::dimensionSet Foam::dimCurrent(0, 0, 0, 0, 0, 1, 0);
const Foam::dimensionSet Foam::dimLuminousIntensity(0, 0, 0, 0, 0, 0, 1);

const Foam::dimensionSet Foam::dimArea(0, 2, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimVolume(0, 3, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimVol(0, 3, 0, 0, 0, 0, 0);

const Foam::dimensionSet Foam::dimDensity(1, -3, 0, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimForce(1, 1, -2, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimEnergy(1, 2, -2, 0, 0, 0, 0);

const Foam::dimensionSet Foam::dimVelocity(0, 1, -1, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimAcceleration(0, 1, -2, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimPressure(1, -1, -2, 0, 0, 0, 0);
const Foam::dimensionSet Foam::dimGasConstant(0, 2, -2, -1, 0, 0, 0);
const Foam::dimensionSet Foam::dimSpecificHeatCapacity(0, 2, -2, -1, 0, 0, 0);


// ************************************************************************* //
