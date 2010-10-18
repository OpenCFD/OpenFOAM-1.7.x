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

#include "LienLeschzinerLowRe.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LienLeschzinerLowRe, 0);
addToRunTimeSelectionTable(RASModel, LienLeschzinerLowRe, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LienLeschzinerLowRe::LienLeschzinerLowRe
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    RASModel(typeName, U, phi, lamTransportModel),

    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.92
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.3
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),
    Am_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Am",
            coeffDict_,
            0.016
        )
    ),
    Aepsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Aepsilon",
            coeffDict_,
            0.263
        )
    ),
    Amu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Amu",
            coeffDict_,
            0.00222
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    epsilon_
    (
        IOobject
        (
            "epsilon",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),

    y_(mesh_),

    yStar_(sqrt(k_)*y_/nu() + SMALL),

    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateLowReNut("nut", mesh_)
    )
{
    nut_ = Cmu_*(scalar(1) - exp(-Am_*yStar_))
       /(scalar(1) - exp(-Aepsilon_*yStar_) + SMALL)*sqr(k_)
       /(epsilon_ + epsilonSmall_);
    nut_.correctBoundaryConditions();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> LienLeschzinerLowRe::R() const
{
    volTensorField gradU = fvc::grad(U_);

    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(gradU),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> LienLeschzinerLowRe::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> LienLeschzinerLowRe::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
    //- (fvc::grad(U) & fvc::grad(nuEff()))
      - fvc::div(nuEff()*fvc::grad(U)().T())
    );
}


bool LienLeschzinerLowRe::read()
{
    if (RASModel::read())
    {
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        Cmu_.readIfPresent(coeffDict());
        kappa_.readIfPresent(coeffDict());
        Am_.readIfPresent(coeffDict());
        Aepsilon_.readIfPresent(coeffDict());
        Amu_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void LienLeschzinerLowRe::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
    }

    scalar Cmu75 = pow(Cmu_.value(), 0.75);

    volTensorField gradU = fvc::grad(U_);

    // generation term
    volScalarField S2 = symm(gradU) && gradU;

    yStar_ = sqrt(k_)*y_/nu() + SMALL;
    volScalarField Rt = sqr(k_)/(nu()*epsilon_);

    volScalarField fMu =
        (scalar(1) - exp(-Am_*yStar_))
       /(scalar(1) - exp(-Aepsilon_*yStar_) + SMALL);

    volScalarField f2 = scalar(1) - 0.3*exp(-sqr(Rt));

    volScalarField G("RASModel::G", Cmu_*fMu*sqr(k_)/epsilon_*S2);


    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        C1_*G*epsilon_/k_
        // E-term
        + C2_*f2*Cmu75*sqrt(k_)
        /(kappa_*y_*(scalar(1) - exp(-Aepsilon_*yStar_)))
       *exp(-Amu_*sqr(yStar_))*epsilon_
      - fvm::Sp(C2_*f2*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

    #include "LienLeschzinerLowReSetWallDissipation.H"
    #include "wallDissipationI.H"

    solve(epsEqn);
    bound(epsilon_, epsilon0_);


    // Turbulent kinetic energy equation

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(), k_)
      ==
        G
      - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, k0_);


    // Re-calculate viscosity
    nut_ = Cmu_*fMu*sqr(k_)/epsilon_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
