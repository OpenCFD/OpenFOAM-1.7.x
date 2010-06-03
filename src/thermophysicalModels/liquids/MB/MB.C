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

\*---------------------------------------------------------------------------*/

#include "MB.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MB, 0);
    addToRunTimeSelectionTable(liquid, MB,);
    addToRunTimeSelectionTable(liquid, MB, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MB::MB()
:
    liquid
    (
        102.133,
        554.5,
        3.4734e+6,
        0.34,
        0.256,
        187.35,
        1.0102e-1,
        375.90,
        5.7373e-30,
        0.3807,
        1.7713e+4
    ),
    rho_(76.6099633, 0.257, 554.5, 0.2772),
    pv_(107.51, -8112.9, -12.77, 9.2919e-06, 2.0),
    hl_(554.5, 508307.794738233, 0.392, 0.0, 0.0, 0.0),
    cp_(1135.77394182096, 2.89818178257762, 0.0, 0.0, 0.0, 0.0),
    h_(-5255966.14542938, 1135.77394182096, 1.44909089128881, 0.0, 0.0, 0.0),
    cpg_(875.329227575808, 2849.22600922327, 1570.0, 2029.70636327142, 678.3),
    B_
    (
        0.00220496803188,
       -2.42184210783978,
       -401045.695318849,
       -2.85079259397061e+17,
       -3.57377145486767e+19
    ),
    mu_(-12.206, 1141.7, 0.15014, 0.0, 0.0),
    mug_(3.733e-07, 0.6177, 256.5, 0.0),
    K_(0.2298, -0.0003002, 0.0, 0.0, 0.0, 0.0),
    Kg_(1333.1, 0.9962, 12317000000.0, 0.0),
    sigma_(554.5, 0.064084, 1.2418, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 102.133, 28.0) // note: Same as nHeptane
{}


Foam::MB::MB
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


Foam::MB::MB(Istream& is)
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
