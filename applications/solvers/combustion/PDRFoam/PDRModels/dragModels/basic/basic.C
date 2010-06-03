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

#include "basic.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDRDragModels
{
    defineTypeNameAndDebug(basic, 0);
    addToRunTimeSelectionTable(PDRDragModel, basic, dictionary);
};
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDRDragModels::basic::basic
(
    const dictionary& PDRProperties,
    const compressible::RASModel& turbulence,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    PDRDragModel(PDRProperties, turbulence, rho, U, phi),
    Csu("Csu", dimless, PDRDragModelCoeffs_.lookup("Csu")),
    Csk("Csk", dimless, PDRDragModelCoeffs_.lookup("Csk")),

    Aw2_
    (
        "Aw2",
        sqr
        (
            volScalarField
            (
                IOobject
                (
                    "Aw",
                    U_.mesh().time().findInstance(polyMesh::meshSubDir, "Aw"),
                    polyMesh::meshSubDir,
                    U_.mesh(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                U_.mesh()
            )
        )
    ),

    CR_
    (
        IOobject
        (
            "CR",
            U_.mesh().time().findInstance(polyMesh::meshSubDir, "CR"),
            polyMesh::meshSubDir,
            U_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh()
    ),

    CT_
    (
        IOobject
        (
            "CT",
            U_.mesh().time().findInstance(polyMesh::meshSubDir, "CT"),
            polyMesh::meshSubDir,
            U_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        U_.mesh()
    )
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::PDRDragModels::basic::~basic()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volSymmTensorField> Foam::PDRDragModels::basic::Dcu() const
{
    const volScalarField& betav = U_.db().lookupObject<volScalarField>("betav");

    return (0.5*rho_)*CR_*mag(U_) + (Csu*I)*betav*turbulence_.muEff()*Aw2_;
}


Foam::tmp<Foam::volScalarField> Foam::PDRDragModels::basic::Gk() const
{
    const volScalarField& betav = U_.db().lookupObject<volScalarField>("betav");

    return
        (0.5*rho_)*mag(U_)*(U_ & CT_ & U_)
      + Csk*betav*turbulence_.muEff()*Aw2_*magSqr(U_);
}


bool Foam::PDRDragModels::basic::read(const dictionary& PDRProperties)
{
    PDRDragModel::read(PDRProperties);

    PDRDragModelCoeffs_.lookup("Csu") >> Csu.value();
    PDRDragModelCoeffs_.lookup("Csk") >> Csk.value();

    return true;
}


// ************************************************************************* //
