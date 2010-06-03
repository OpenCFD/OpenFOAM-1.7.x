/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#ifndef makeReactingParcelPhaseChangeModels_H
#define makeReactingParcelPhaseChangeModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"
#include "ReactingCloud.H"

#include "NoPhaseChange.H"
#include "LiquidEvaporation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactingPhaseChangeModels(ParcelType)                             \
                                                                              \
    makeReactingPhaseChangeModelThermoType                                    \
    (                                                                         \
        ParcelType,                                                           \
        constGasThermoPhysics                                                 \
    );                                                                        \
                                                                              \
    makeReactingPhaseChangeModelThermoType                                    \
    (                                                                         \
        ParcelType,                                                           \
        gasThermoPhysics                                                      \
    );                                                                        \
                                                                              \
    makeReactingPhaseChangeModelThermoType                                    \
    (                                                                         \
        ParcelType,                                                           \
        icoPoly8ThermoPhysics                                                 \
    );


#define makeReactingPhaseChangeModelThermoType(ParcelType, ThermoType)        \
                                                                              \
    makePhaseChangeModel(ReactingCloud<ParcelType<ThermoType> >);             \
                                                                              \
    makePhaseChangeModelThermoType                                            \
    (                                                                         \
        NoPhaseChange,                                                        \
        ReactingCloud,                                                        \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makePhaseChangeModelThermoType                                            \
    (                                                                         \
        LiquidEvaporation,                                                    \
        ReactingCloud,                                                        \
        ParcelType,                                                           \
        ThermoType                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
