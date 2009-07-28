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

#include "C2H6.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C2H6, 0);
    addToRunTimeSelectionTable(liquid, C2H6,);
    addToRunTimeSelectionTable(liquid, C2H6, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C2H6::C2H6()
:
    liquid
    (
        30.070,
        305.32,
        4.872e+6,
        0.14550,
        0.279,
        90.35,
        1.13,
        184.55,
        0.0,
        0.0995,
        1.24e+4
    ),
    rho_(57.499854, 0.27937, 305.32, 0.29187),
    pv_(51.857, -2598.7, -5.1283, 1.4913e-05, 2.0),
    hl_(305.32, 701396.740937812, 0.60646, -0.55492, 0.32799, 0.0),
    cp_
    (
        305.32,
        8.02554965861611,
        2983.63817758563,
        167.548325566287,
       -343.93389207094
    ),
    h_(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    cpg_(1341.07083471899, 4463.58496840705, 1655.5, 2435.08480212837, 752.87),
    B_
    (
        0.00269205187894912,
       -2.05221150648487,
       -47721.9820419022,
        2.24808779514466e+15,
       -3.23910874625873e+17
    ),
    mu_(-3.4134, 197.05, -1.2193, -9.2023e-26, 10.0),
    mug_(2.5906e-07, 0.67988, 98.902, 0.0),
    K_(0.35758, -0.0011458, 6.1866e-07, 0.0, 0.0, 0.0),
    Kg_(7.3869e-05, 1.1689, 500.73, 0.0),
    sigma_(305.32, 0.048643, 1.1981, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 30.070, 28) // note: Same as nHeptane
{}


Foam::C2H6::C2H6
(
    const liquid& l,
    const NSRDSfunc5& density,
    const NSRDSfunc1& vapourPressure,
    const NSRDSfunc6& heatOfVapourisation,
    const NSRDSfunc14& heatCapacity,
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


Foam::C2H6::C2H6(Istream& is)
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
