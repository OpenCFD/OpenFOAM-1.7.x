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

#include "IC8H18.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(IC8H18, 0);
    addToRunTimeSelectionTable(liquid, IC8H18,);
    addToRunTimeSelectionTable(liquid, IC8H18, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IC8H18::IC8H18()
:
    liquid
    (
        114.231,
        543.96,
        2.5676e+6,
        0.468,
        0.266,
        165.78,
        1.4464e-2,
        372.39,
        0.0,
        0.3031,
        1.4051e+4
    ),
    rho_(67.2363666, 0.27373, 543.96, 0.2846),
    pv_(120.81, -7550, -16.111, 0.017099, 1.0),
    hl_(543.96, 375379.713037617, 0.1549, 0.138, 0.0666, 0.0),
    cp_
    (
        1219.89652546156,
        1.67205049417409,
        0.00414073237562483,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
       -2743275.10767575,
        1219.89652546156,
        0.836025247087043,
        0.00138024412520828,
        0.0,
        0.0
    ),
    cpg_(997.10236275617, 4627.4653990598, 1594, 2933.52942721328, 677.94),
    B_
    (
        0.00234936225718063,
       -2.83381919093766,
       -413154.047500241,
       -3.49703670632315e+17,
        3.13750207912038e+19
    ),
    mu_(-15.811, 1282.5, 0.67791, -3.8617e-28, 10.0),
    mug_(1.107e-07, 0.746, 72.4, 0.0),
    K_(0.1508, -0.0001712, 0.0, 0.0, 0.0, 0.0),
    Kg_(1.758e-05, 1.3114, 392.9, 0.0),
    sigma_(543.96, 0.047434, 1.1975, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 114.231, 28.0) // note: Same as nHeptane
{}


Foam::IC8H18::IC8H18
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


Foam::IC8H18::IC8H18(Istream& is)
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
