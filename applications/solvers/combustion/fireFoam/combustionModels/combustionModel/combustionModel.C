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

#include "combustionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(combustionModel, 0);
    defineRunTimeSelectionTable(combustionModel, dictionary);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModel::combustionModel
(
    const dictionary& combustionProperties,
    const hsCombustionThermo& thermo,
    const compressible::turbulenceModel& turbulence,
    const surfaceScalarField& phi,
    const volScalarField& rho
)
:
    combustionModelCoeffs_
    (
        combustionProperties.subDict
        (
            word(combustionProperties.lookup("combustionModel")) + "Coeffs"
        )
    ),
    thermo_(thermo),
    turbulence_(turbulence),
    mesh_(phi.mesh()),
    phi_(phi),
    rho_(rho),
    stoicRatio_(thermo.lookup("stoichiometricAirFuelMassRatio")),
    s_(thermo.lookup("stoichiometricOxygenFuelMassRatio")),
    qFuel_(thermo_.lookup("qFuel")),
    composition_(thermo.composition())
{}



Foam::tmp<Foam::volScalarField> Foam::combustionModel::combustionModel::Sh
(
    const volScalarField& Rfu
) const
{
    return Rfu*qFuel_;
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::combustionModel::~combustionModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::combustionModel::read(const dictionary& combustionProperties)
{
    combustionModelCoeffs_ = combustionProperties.subDict(type() + "Coeffs");

    return true;
}


// ************************************************************************* //
