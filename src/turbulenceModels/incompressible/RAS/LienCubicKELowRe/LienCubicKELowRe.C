/*---------------------------------------------------------------------------* \
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

#include "LienCubicKELowRe.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(LienCubicKELowRe, 0);
addToRunTimeSelectionTable(RASModel, LienCubicKELowRe, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

LienCubicKELowRe::LienCubicKELowRe
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
    A1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A1",
            coeffDict_,
            1.25
        )
    ),
    A2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A2",
            coeffDict_,
            1000.0
        )
    ),
    Ctau1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ctau1",
            coeffDict_,
            -4.0
        )
    ),
    Ctau2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ctau2",
            coeffDict_,
            13.0
        )
    ),
    Ctau3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ctau3",
            coeffDict_,
            -2.0
        )
    ),
    alphaKsi_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaKsi",
            coeffDict_,
            0.9
        )
    ),
    CmuWall_
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

    gradU_(fvc::grad(U)),
    eta_(k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU_ + gradU_.T())))),
    ksi_(k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU_ - gradU_.T())))),
    Cmu_(2.0/(3.0*(A1_ + eta_ + alphaKsi_*ksi_))),
    fEta_(A2_ + pow(eta_, 3.0)),

    C5viscosity_
    (
        -2.0*pow(Cmu_, 3.0)*pow(k_, 4.0)/pow(epsilon_, 3.0)
       *(magSqr(gradU_ + gradU_.T()) - magSqr(gradU_ - gradU_.T()))
    ),

    yStar_(sqrt(k_)*y_/nu() + SMALL),

    nut_
    (
        Cmu_
       *(
            scalar(1) - exp(-Am_*yStar_))
           /(scalar(1) - exp(-Aepsilon_*yStar_) + SMALL
        )
       *sqr(k_)/(epsilon_ + epsilonSmall_)
        // cubic term C5, implicit part
      + max
        (
            C5viscosity_,
            dimensionedScalar("0", C5viscosity_.dimensions(), 0.0)
        )
    ),
    // turbulent viscosity, with implicit part of C5

    nonlinearStress_
    (
        "nonlinearStress",
        symm
        (
            // quadratic terms
            pow(k_, 3.0)/sqr(epsilon_)
           *(
                Ctau1_/fEta_
               *(
                    (gradU_ & gradU_)
                  + (gradU_ & gradU_)().T()
                )
              + Ctau2_/fEta_*(gradU_ & gradU_.T())
              + Ctau3_/fEta_*(gradU_.T() & gradU_)
            )
            // cubic term C4
          - 20.0*pow(k_, 4.0)/pow(epsilon_, 3.0)
           *pow(Cmu_, 3.0)
           *(
                ((gradU_ & gradU_) & gradU_.T())
              + ((gradU_ & gradU_.T()) & gradU_.T())
              - ((gradU_.T() & gradU_) & gradU_)
              - ((gradU_.T() & gradU_.T()) & gradU_)
            )
            // cubic term C5, explicit part
          + min
            (
                C5viscosity_,
                dimensionedScalar("0", C5viscosity_.dimensions(), 0.0)
            )*gradU_
        )
    )
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> LienCubicKELowRe::R() const
{
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
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(gradU_) + nonlinearStress_,
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> LienCubicKELowRe::devReff() const
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
           -nuEff()*dev(twoSymm(fvc::grad(U_))) + nonlinearStress_
        )
    );
}


tmp<fvVectorMatrix> LienCubicKELowRe::divDevReff(volVectorField& U) const
{
    return
    (
        fvc::div(nonlinearStress_)
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool LienCubicKELowRe::read()
{
    if (RASModel::read())
    {
        C1_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        A1_.readIfPresent(coeffDict());
        A2_.readIfPresent(coeffDict());
        Ctau1_.readIfPresent(coeffDict());
        Ctau2_.readIfPresent(coeffDict());
        Ctau3_.readIfPresent(coeffDict());
        alphaKsi_.readIfPresent(coeffDict());
        CmuWall_.readIfPresent(coeffDict());
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


void LienCubicKELowRe::correct()
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

    gradU_ = fvc::grad(U_);

    // generation term
    volScalarField S2 = symm(gradU_) && gradU_;

    yStar_ = sqrt(k_)*y_/nu() + SMALL;
    volScalarField Rt = sqr(k_)/(nu()*epsilon_);

    volScalarField fMu =
        (scalar(1) - exp(-Am_*yStar_))
       /(scalar(1) - exp(-Aepsilon_*yStar_) + SMALL);

    volScalarField f2 = scalar(1) - 0.3*exp(-sqr(Rt));

    volScalarField G =
        Cmu_*fMu*sqr(k_)/epsilon_*S2 - (nonlinearStress_ && gradU_);

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
      ==
        C1_*G*epsilon_/k_
        // E-term
      + C2_*f2*pow(Cmu_, 0.75)*pow(k_, scalar(0.5))
       /(kappa_*y_*(scalar(1) - exp(-Aepsilon_*yStar_)))
       *exp(-Amu_*sqr(yStar_))*epsilon_
      - fvm::Sp(C2_*f2*epsilon_/k_, epsilon_)
    );

    epsEqn().relax();

#   include "LienCubicKELowReSetWallDissipation.H"
#   include "wallDissipationI.H"

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

    eta_ = k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU_ + gradU_.T())));
    ksi_ = k_/epsilon_*sqrt(2.0*magSqr(0.5*(gradU_ - gradU_.T())));
    Cmu_ = 2.0/(3.0*(A1_ + eta_ + alphaKsi_*ksi_));
    fEta_ = A2_ + pow(eta_, 3.0);

    C5viscosity_ =
      - 2.0*pow(Cmu_, 3.0)*pow(k_, 4.0)/pow(epsilon_, 3.0)
       *(magSqr(gradU_ + gradU_.T()) - magSqr(gradU_ - gradU_.T()));

    nut_ =
        Cmu_*fMu*sqr(k_)/epsilon_
        // C5 term, implicit
      + max
        (
            C5viscosity_,
            dimensionedScalar("0", C5viscosity_.dimensions(), 0.0)
        );

    nonlinearStress_ = symm
    (
        // quadratic terms
        pow(k_, 3.0)/sqr(epsilon_)
       *(
            Ctau1_/fEta_
           *(
                (gradU_ & gradU_)
              + (gradU_ & gradU_)().T()
            )
          + Ctau2_/fEta_*(gradU_ & gradU_.T())
          + Ctau3_/fEta_*(gradU_.T() & gradU_)
        )
        // cubic term C4
      - 20.0*pow(k_, 4.0)/pow(epsilon_, 3.0)
       *pow(Cmu_, 3.0)
       *(
            ((gradU_ & gradU_) & gradU_.T())
          + ((gradU_ & gradU_.T()) & gradU_.T())
          - ((gradU_.T() & gradU_) & gradU_)
          - ((gradU_.T() & gradU_.T()) & gradU_)
        )
        // cubic term C5, explicit part
      + min
        (
            C5viscosity_,
            dimensionedScalar("0", C5viscosity_.dimensions(), 0.0)
        )*gradU_
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
