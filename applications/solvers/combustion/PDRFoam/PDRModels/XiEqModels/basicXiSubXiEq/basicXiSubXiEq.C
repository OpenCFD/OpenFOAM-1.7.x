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

#include "basicXiSubXiEq.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace XiEqModels
{
    defineTypeNameAndDebug(basicSubGrid, 0);
    addToRunTimeSelectionTable(XiEqModel, basicSubGrid, dictionary);
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XiEqModels::basicSubGrid::basicSubGrid
(
    const dictionary& XiEqProperties,
    const hhuCombustionThermo& thermo,
    const compressible::RASModel& turbulence,
    const volScalarField& Su
)
:
    XiEqModel(XiEqProperties, thermo, turbulence, Su),

    N_
    (
        IOobject
        (
            "N",
            Su.mesh().time().findInstance(polyMesh::meshSubDir, "N"),
            polyMesh::meshSubDir,
            Su.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        Su.mesh()
    ),

    ns_
    (
        IOobject
        (
            "ns",
            Su.mesh().time().findInstance(polyMesh::meshSubDir, "ns"),
            polyMesh::meshSubDir,
            Su.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        Su.mesh()
    ),

    B_
    (
        IOobject
        (
            "B",
            Su.mesh().time().findInstance(polyMesh::meshSubDir, "B"),
            polyMesh::meshSubDir,
            Su.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        Su.mesh()
    ),

    Lobs_
    (
        IOobject
        (
            "Lobs",
            Su.mesh().time().findInstance(polyMesh::meshSubDir, "Lobs"),
            polyMesh::meshSubDir,
            Su.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        Su.mesh()
    ),

    XiEqModel_(XiEqModel::New(XiEqModelCoeffs_, thermo, turbulence, Su))
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::XiEqModels::basicSubGrid::~basicSubGrid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::XiEqModels::basicSubGrid::XiEq() const
{
    const objectRegistry& db = Su_.db();
    const volVectorField& U = db.lookupObject<volVectorField>("U");

    volScalarField magU = mag(U);
    volVectorField Uhat =
    U/(mag(U) + dimensionedScalar("Usmall", U.dimensions(), 1e-4));

    volScalarField n = max(N_ - (Uhat & ns_ & Uhat), scalar(1e-4));

    volScalarField b = (Uhat & B_ & Uhat)/n;

    volScalarField up = sqrt((2.0/3.0)*turbulence_.k());

    volScalarField XiSubEq =
        scalar(1)
      + max(2.2*sqrt(b), min(0.34*magU/up, scalar(1.6)))
       *min(0.25*n, scalar(1));

    return XiSubEq*XiEqModel_->XiEq();
}


bool Foam::XiEqModels::basicSubGrid::read(const dictionary& XiEqProperties)
{
    XiEqModel::read(XiEqProperties);

    return XiEqModel_->read(XiEqModelCoeffs_);
}


// ************************************************************************* //
