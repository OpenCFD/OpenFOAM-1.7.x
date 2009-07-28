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

#include "bC10H7CH3.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(bC10H7CH3, 0);
    addToRunTimeSelectionTable(liquid, bC10H7CH3,);
    addToRunTimeSelectionTable(liquid, bC10H7CH3, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::bC10H7CH3::bC10H7CH3()
:
    liquid
    (
        142.2,
        761.0,
        3.25e+6,
        0.507,
        0.260,
        307.73,
        1.7374e+1,
        514.20,
        1.4010e-30,
        0.3459,
        1.987e+4
    ),
    rho_(67.36014, 0.23843, 761, 0.2559),
    pv_(134.31, -12103, -16.195, 6.9659e-06, 2),
    hl_(761.0, 513150.492264416, 0.4044, 0.0, 0.0, 0.0),
    cp_(811.322081575246, 2.30225035161744, 0.0008628691983122, 0.0, 0.0, 0.0),
    h_
    (
        45001.2311880177,
        811.322081575246,
        1.15112517580872,
        0.000287623066104079,
        0.0,
        0.0
    ),
    cpg_(760.126582278481, 2699.08579465542, 1564.1, 1994.51476793249, 727.49),
    B_
    (
        0.00229430379746835,
       -3.53720112517581,
       -1067158.93108298,
        2.29746835443038e+18,
       -2.68438818565401e+21
    ),
    mu_(-63.276, 4219, 7.5549, 0.0, 0.0),
    mug_(2.1791e-06, 0.3717, 712.53, 0.0),
    K_(0.1962, -0.00018414, 0.0, 0.0, 0.0, 0.0),
    Kg_(0.4477, -0.1282, -345.89, 2340100),
    sigma_(761.0, 0.066442, 1.2634, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 142.2, 28) // note: Same as nHeptane
{}


Foam::bC10H7CH3::bC10H7CH3
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


Foam::bC10H7CH3::bC10H7CH3(Istream& is)
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
