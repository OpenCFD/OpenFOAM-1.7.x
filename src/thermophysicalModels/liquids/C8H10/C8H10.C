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

#include "C8H10.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(C8H10, 0);
    addToRunTimeSelectionTable(liquid, C8H10,);
    addToRunTimeSelectionTable(liquid, C8H10, Istream);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::C8H10::C8H10()
:
    liquid
    (
        106.167,
        617.17,
        3.6094e+6,
        0.37381,
        0.263,
        178.15,
        4.038e-3,
        409.35,
        1.9680e-30,
        0.3036,
        1.8043e+4
    ),
    rho_(76.3765398, 0.26438, 617.17, 0.2921),
    pv_(88.246, -7691.1, -9.797, 5.931e-06, 2.0),
    hl_(617.17, 516167.924119547, 0.3882, 0.0, 0.0, 0.0),
    cp_
    (
        818.521762883005,
        6.66873887366131,
       -0.0248005500767658,
        4.23860521631015e-05,
        0.0,
        0.0
    ),
    h_
    (
       -524002.612929508,
        818.521762883005,
        3.33436943683065,
       -0.00826685002558862,
        1.05965130407754e-05,
        0.0
    ),
    cpg_(738.835984816374, 3201.5598067196, 1559, 2285.07916772632, -702.0),
    B_
    (
        0.00165776559571242,
       -2.77958310962917,
       -388067.855359952,
       -5.86905535618412e+18,
        1.58052878954854e+21
    ),
    mu_(-10.452, 1048.4, -0.0715, 0.0, 0.0),
    mug_(1.2e-06, 0.4518, 439.0, 0.0),
    K_(0.20149, -0.00023988, 0.0, 0.0, 0.0, 0.0),
    Kg_(1.708e-05, 1.319, 565.6, 0.0),
    sigma_(617.17, 0.066, 1.268, 0.0, 0.0, 0.0),
    D_(147.18, 20.1, 106.167, 28.0) // note: Same as nHeptane
{}


Foam::C8H10::C8H10
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


Foam::C8H10::C8H10(Istream& is)
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
