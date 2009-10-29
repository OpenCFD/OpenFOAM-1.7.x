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

#include "radiationModel.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(radiationModel, 0);
        defineRunTimeSelectionTable(radiationModel, dictionary);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::radiationModel::radiationModel(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    radiation_(false),
    coeffs_(dictionary::null),
    solverFreq_(0),
    absorptionEmission_(NULL),
    scatter_(NULL)
{}


Foam::radiation::radiationModel::radiationModel
(
    const word& type,
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    radiation_(lookup("radiation")),
    coeffs_(subDict(type + "Coeffs")),
    solverFreq_(readLabel(lookup("solverFreq"))),
    absorptionEmission_(absorptionEmissionModel::New(*this, mesh_)),
    scatter_(scatterModel::New(*this, mesh_))
{
    solverFreq_ = max(1, solverFreq_);
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::radiationModel::~radiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::radiationModel::read()
{
    if (regIOobject::read())
    {
        lookup("radiation") >> radiation_;
        coeffs_ = subDict(type() + "Coeffs");

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::radiation::radiationModel::correct()
{
    if (!radiation_)
    {
        return;
    }

    if (time_.timeIndex() % solverFreq_ == 0)
    {
        calculate();
    }
}


Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModel::Sh
(
    basicThermo& thermo
) const
{
    volScalarField& h = thermo.h();
    const volScalarField cp = thermo.Cp();
    const volScalarField T3 = pow3(T_);

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/cp, h)
      - Rp()*T3*(T_ - 4.0*h/cp)
    );
}

Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModel::Shs
(
    basicThermo& thermo
) const
{
    volScalarField& hs = thermo.hs();
    const volScalarField cp = thermo.Cp();
    const volScalarField T3 = pow3(T_);

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/cp, hs)
      - Rp()*T3*(T_ - 4.0*hs/cp)
    );
}

// ************************************************************************* //
