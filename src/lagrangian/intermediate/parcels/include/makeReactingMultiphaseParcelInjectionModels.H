/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#ifndef makeReactingMultiphaseParcelInjectionModels_H
#define makeReactingMultiphaseParcelInjectionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"
#include "KinematicCloud.H"

#include "ConeInjection.H"
#include "ConeInjectionMP.H"
#include "FieldActivatedInjection.H"
#include "ManualInjection.H"
#include "NoInjection.H"
#include "PatchInjection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactingMultiphaseInjectionModels(ParcelType)                     \
                                                                              \
    makeReactingMultiphaseInjectionModelThermoType                            \
    (                                                                         \
        ParcelType,                                                           \
        constGasThermoPhysics                                                 \
    );                                                                        \
                                                                              \
    makeReactingMultiphaseInjectionModelThermoType                            \
    (                                                                         \
        ParcelType,                                                           \
        gasThermoPhysics                                                      \
    );                                                                        \
                                                                              \
    makeReactingMultiphaseInjectionModelThermoType                            \
    (                                                                         \
        ParcelType,                                                           \
        icoPoly8ThermoPhysics                                                 \
    );


#define makeReactingMultiphaseInjectionModelThermoType(ParcelType, ThermoType)\
                                                                              \
    makeInjectionModel(KinematicCloud<ParcelType<ThermoType> >);              \
                                                                              \
    makeInjectionModelThermoType                                              \
    (                                                                         \
        ConeInjection,                                                        \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeInjectionModelThermoType                                              \
    (                                                                         \
        ConeInjectionMP,                                                      \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeInjectionModelThermoType                                              \
    (                                                                         \
        FieldActivatedInjection,                                              \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeInjectionModelThermoType                                              \
    (                                                                         \
        ManualInjection,                                                      \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeInjectionModelThermoType                                              \
    (                                                                         \
        NoInjection,                                                          \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeInjectionModelThermoType                                              \
    (                                                                         \
        PatchInjection,                                                       \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
