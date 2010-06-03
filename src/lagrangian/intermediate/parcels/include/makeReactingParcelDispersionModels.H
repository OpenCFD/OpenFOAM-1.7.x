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

#ifndef makeReactingParcelDispersionModels_H
#define makeReactingParcelDispersionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"
#include "KinematicCloud.H"

#include "NoDispersion.H"
#include "GradientDispersionRAS.H"
#include "StochasticDispersionRAS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactingDispersionModels(ParcelType)                              \
                                                                              \
    makeReactingDispersionModelThermoType                                     \
    (                                                                         \
        ParcelType,                                                           \
        constGasThermoPhysics                                                 \
    );                                                                        \
                                                                              \
    makeReactingDispersionModelThermoType                                     \
    (                                                                         \
        ParcelType,                                                           \
        gasThermoPhysics                                                      \
    );                                                                        \
                                                                              \
    makeReactingDispersionModelThermoType                                     \
    (                                                                         \
        ParcelType,                                                           \
        icoPoly8ThermoPhysics                                                 \
    );


#define makeReactingDispersionModelThermoType(ParcelType, ThermoType)         \
                                                                              \
    makeDispersionModel(KinematicCloud<ParcelType<ThermoType> >);             \
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        DispersionRASModel<KinematicCloud<ParcelType<ThermoType> > >,         \
        0                                                                     \
    );                                                                        \
                                                                              \
    makeDispersionModelThermoType                                             \
    (                                                                         \
        NoDispersion,                                                         \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeDispersionModelThermoType                                             \
    (                                                                         \
        GradientDispersionRAS,                                                \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeDispersionModelThermoType                                             \
    (                                                                         \
        StochasticDispersionRAS,                                              \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
