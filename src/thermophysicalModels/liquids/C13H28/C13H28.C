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

#include "C13H28.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C13H28, 0);
    addToRunTimeSelectionTable(liquid, C13H28,);
    addToRunTimeSelectionTable(liquid, C13H28, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C13H28::C13H28()
:
    liquid
    (
        184.365,
        675.80,
        1.7225e+6,
        0.77,
        0.236,
        267.76,
        3.801e-1,
        508.62,
        0.0,
        0.6186,
        1.5901e+4
    ),
    rho_(59.513022, 0.2504, 675.8, 0.312),
    pv_(118.27, -11432, -13.769, 5.9641e-06, 2.0),
    hl_(675.80, 444227.48352453, 0.4162, 0.0, 0.0, 0.0),
    cp_
    (
        4275.05220622135,
       -16.6539202126217,
        0.0325755973205326,
        0.0,
        0.0,
        0.0
    ),
    h_
    (
       -2860442.0545124,
        4275.05220622135,
       -8.32696010631085,
        0.0108585324401775,
        0.0,
        0.0
    ),
    cpg_(1136.87522035093, 3641.14663846175, -1443, 2277.00485450058, -683.0),
    B_
    (
        0.00246321156401703,
       -2.66601578390692,
       -1249532.17801643,
       -1.0460770753668e+19,
        1.90117430097904e+21
    ),
    mu_(-23.341, 2121.9, 1.7208, 0.0, 0.0),
    mug_(3.5585e-08, 0.8987, 165.3, 0.0),
    K_(0.1981, -0.0002046, 0.0, 0.0, 0.0, 0.0),
    Kg_(5.3701e-06, 1.4751, 599.09, 0.0),
    sigma_(675.80, 0.05561, 1.3361, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 184.365, 28.0) // note: Same as nHeptane
{}


Foam::C13H28::C13H28
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


Foam::C13H28::C13H28(Istream& is)
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
