/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "instabilityG.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiGModels
{
    defineTypeNameAndDebug(instabilityG, 0);
    addToRunTimeSelectionTable(XiGModel, instabilityG, dictionary);
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiGModels::instabilityG::instabilityG
(
    const dictionary& XiGProperties,
    const hhuCombustionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiGModel(XiGProperties, thermo, turbulence, Su),
    GIn(XiGModelCoeffs_.lookup("GIn")),
    lambdaIn(XiGModelCoeffs_.lookup("lambdaIn")),
    XiGModel_(XiGModel::New(XiGModelCoeffs_, thermo, turbulence, Su))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::XiGModels::instabilityG::~instabilityG()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiGModels::instabilityG::G() const
{
    volScalarField turbXiG = XiGModel_->G();
    return GIn*GIn/(GIn + turbXiG) + turbXiG;
}


Foam::tmp<Foam::volScalarField> Foam::XiGModels::instabilityG::Db() const
{
    const objectRegistry& db = Su_.db();
    const volScalarField& Xi = db.lookupObject<volScalarField>("Xi");
    const volScalarField& rho = db.lookupObject<volScalarField>("rho");
    const volScalarField& mgb = db.lookupObject<volScalarField>("mgb");

    return XiGModel_->Db()
        + rho*Su_*(Xi - 1.0)*mgb*(0.5*lambdaIn)/(mgb + 1.0/lambdaIn);
}


bool Foam::XiGModels::instabilityG::read(const dictionary& XiGProperties)
{
    XiGModel::read(XiGProperties);

    XiGModelCoeffs_.lookup("GIn") >> GIn;
    XiGModelCoeffs_.lookup("lambdaIn") >> lambdaIn;

    return true;
}


// ************************************************************************* //
