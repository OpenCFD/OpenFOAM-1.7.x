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

#ifndef makeReactingMultiphaseParcelSurfaceReactionModels_H
#define makeReactingMultiphaseParcelSurfaceReactionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"
#include "ReactingMultiphaseCloud.H"

#include "NoSurfaceReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactingMultiphaseSurfaceReactionModels(ParcelType)               \
                                                                              \
    makeReactingMultiphaseSurfaceReactionModelThermoType                      \
    (                                                                         \
        ParcelType,                                                           \
        constGasThermoPhysics                                                 \
    );                                                                        \
    makeReactingMultiphaseSurfaceReactionModelThermoType                      \
    (                                                                         \
        ParcelType,                                                           \
        gasThermoPhysics                                                      \
    );                                                                        \
    makeReactingMultiphaseSurfaceReactionModelThermoType                      \
    (                                                                         \
        ParcelType,                                                           \
        icoPoly8ThermoPhysics                                                 \
    );


#define makeReactingMultiphaseSurfaceReactionModelThermoType(ParcelType, ThermoType)\
                                                                              \
    makeSurfaceReactionModel                                                  \
    (                                                                         \
        ReactingMultiphaseCloud<ParcelType<ThermoType> >                      \
    );                                                                        \
                                                                              \
    makeSurfaceReactionModelThermoType                                        \
    (                                                                         \
        NoSurfaceReaction,                                                    \
        ReactingMultiphaseCloud,                                              \
        ParcelType,                                                           \
        ThermoType                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
