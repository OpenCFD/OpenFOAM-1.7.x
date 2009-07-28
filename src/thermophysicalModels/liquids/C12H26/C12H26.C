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

#include "C12H26.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C12H26, 0);
    addToRunTimeSelectionTable(liquid, C12H26,);
    addToRunTimeSelectionTable(liquid, C12H26, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C12H26::C12H26()
:
    liquid
    (
        170.338,
        658.0,
        1.82e+6,
        0.716,
        0.238,
        263.57,
        6.152e-1,
        489.47,
        0.0,
        0.5764,
        1.59e+4
    ),
    rho_(60.53982858, 0.25511, 658.0, 0.29368),
    pv_(137.47, -11976.0, -16.698, 8.0906e-06, 2.0),
    hl_(658.0, 454020.829174935, 0.40681, 0.0, 0.0, 0.0),
    cp_(2983.53861146661, -8.0352006011577, 0.018207916025784, 0.0, 0.0, 0.0),
    h_
    (
       -2755166.83820769,
        2983.53861146661,
       -4.01760030057885,
        0.00606930534192801,
        0.0,
        0.0
    ),
    cpg_(1250.16144371778, 3894.02247296552, 1715.5, 2650.67101879792, 777.5),
    B_
    (
        0.00516619896910848,
       -6.40491258556517,
       -295295.236529723,
       -3.22147729807794e+19,
        8.78195117941974e+21
    ),
    mu_(-20.607, 1943, 1.3205, 0.0, 0.0),
    mug_(6.344e-08, 0.8287, 219.5, 0.0),
    K_(0.2047, -0.0002326, 0.0, 0.0, 0.0, 0.0),
    Kg_(5.719e-06, 1.4699, 579.4, 0.0),
    sigma_(658.0, 0.055493, 1.3262, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 170.338, 28.0) // note: Same as nHeptane
{}


Foam::C12H26::C12H26
(
    const liquid& l,
    const NSRDSfunc5& density,
    const NSRDSfunc1& vapourPressure,
    const NSRDSfunc6& heatOfVapourisation,
    const NSRDSfunc0& heatCapacity,
    const NSRDSfunc0& enthalpy,
    const NSRDSfunc7& idealGasHeatCapacity,
    const NSRDSfunc4& secondVirialCoeff,
    const NSRDSfunc1& dynamicViscosity,
    const NSRDSfunc2& vapourDynamicViscosity,
    const NSRDSfunc0& thermalConductivity,
    const NSRDSfunc2& vapourThermalConductivity,
    const NSRDSfunc6& surfaceTension,
    const APIdiffCoefFunc& vapourDiffussivity
)
:
    liquid(l),
    rho_(density),
    pv_(vapourPressure),
    hl_(heatOfVapourisation),
    cp_(heatCapacity),
    h_(enthalpy),
    cpg_(idealGasHeatCapacity),
    B_(secondVirialCoeff),
    mu_(dynamicViscosity),
    mug_(vapourDynamicViscosity),
    K_(thermalConductivity),
    Kg_(vapourThermalConductivity),
    sigma_(surfaceTension),
    D_(vapourDiffussivity)
{}


Foam::C12H26::C12H26(Istream& is)
:
    liquid(is),
    rho_(is),
    pv_(is),
    hl_(is),
    cp_(is),
    h_(is),
    cpg_(is),
    B_(is),
    mu_(is),
    mug_(is),
    K_(is),
    Kg_(is),
    sigma_(is),
    D_(is)
{}


// ************************************************************************* //
