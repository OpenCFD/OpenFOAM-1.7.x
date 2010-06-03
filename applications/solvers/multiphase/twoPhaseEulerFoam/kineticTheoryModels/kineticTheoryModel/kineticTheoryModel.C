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

#include "kineticTheoryModel.H"
#include "surfaceInterpolate.H"
#include "mathematicalConstants.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::kineticTheoryModel
(
    const Foam::phaseModel& phasea,
    const Foam::volVectorField& Ub,
    const Foam::volScalarField& alpha,
    const Foam::dragModel& draga
)
:
    phasea_(phasea),
    Ua_(phasea.U()),
    Ub_(Ub),
    alpha_(alpha),
    phia_(phasea.phi()),
    draga_(draga),

    rhoa_(phasea.rho()),
    da_(phasea.d()),
    nua_(phasea.nu()),

    kineticTheoryProperties_
    (
        IOobject
        (
            "kineticTheoryProperties",
            Ua_.time().constant(),
            Ua_.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    kineticTheory_(kineticTheoryProperties_.lookup("kineticTheory")),
    equilibrium_(kineticTheoryProperties_.lookup("equilibrium")),
    viscosityModel_
    (
        kineticTheoryModels::viscosityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    conductivityModel_
    (
        conductivityModel::New
        (
            kineticTheoryProperties_
        )
    ),
    radialModel_
    (
        radialModel::New
        (
            kineticTheoryProperties_
        )
    ),
    granularPressureModel_
    (
        granularPressureModel::New
        (
            kineticTheoryProperties_
        )
    ),
    frictionalStressModel_
    (
        frictionalStressModel::New
        (
            kineticTheoryProperties_
        )
    ),
    e_(kineticTheoryProperties_.lookup("e")),
    alphaMax_(kineticTheoryProperties_.lookup("alphaMax")),
    alphaMinFriction_(kineticTheoryProperties_.lookup("alphaMinFriction")),
    Fr_(kineticTheoryProperties_.lookup("Fr")),
    eta_(kineticTheoryProperties_.lookup("eta")),
    p_(kineticTheoryProperties_.lookup("p")),
    phi_(dimensionedScalar(kineticTheoryProperties_.lookup("phi"))*M_PI/180.0),
    Theta_
    (
        IOobject
        (
            "Theta",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        Ua_.mesh()
    ),
    mua_
    (
        IOobject
        (
            "mua",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    pa_
    (
        IOobject
        (
            "pa",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -2, 0, 0), 0.0)
    ),
    kappa_
    (
        IOobject
        (
            "kappa",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(1, -1, -1, 0, 0), 0.0)
    ),
    gs0_
    (
        IOobject
        (
            "gs0",
            Ua_.time().timeName(),
            Ua_.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Ua_.mesh(),
        dimensionedScalar("zero", dimensionSet(0, 0, 0, 0, 0), 1.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::kineticTheoryModel::~kineticTheoryModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kineticTheoryModel::solve()
{
    if (!kineticTheory_)
    {
        return;
    }

    word scheme("div(phi,Theta)");

    volScalarField alpha = alpha_;
    alpha.max(1.0e-6);
    const scalar sqrtPi = sqrt(mathematicalConstant::pi);

    surfaceScalarField phi = 1.5*rhoa_*phia_*fvc::interpolate(alpha_);

    volTensorField dU = fvc::grad(Ua_);
    volTensorField dUT = dU.T();
    volTensorField D = 0.5*(dU + dUT);

    // NB, drag = K*alpha*beta,
    // (the alpha and beta has been extracted from the drag function for
    // numerical reasons)
    volScalarField Ur = mag(Ua_ - Ub_);
    volScalarField betaPrim = alpha_*(1.0 - alpha_)*draga_.K(Ur);

    // Calculating the radial distribution function (solid volume fraction is
    //  limited close to the packing limit, but this needs improvements)
    //  The solution is higly unstable close to the packing limit.
    gs0_ = radialModel_->g0(min(alpha, alphaMax_-1.0e-2), alphaMax_);

    // particle pressure - coefficient in front of Theta (Eq. 3.22, p. 45)
    volScalarField PsCoeff =
        granularPressureModel_->granularPressureCoeff(alpha_,gs0_,rhoa_,e_ );

    dimensionedScalar Tsmall
    (
        "small",
        dimensionSet(0,2,-2,0,0,0,0),
        1.0e-6
    );

    dimensionedScalar TsmallSqrt = sqrt(Tsmall);
    volScalarField ThetaSqrt = sqrt(Theta_);

    // 'thermal' conductivity (Table 3.3, p. 49)
    kappa_ = conductivityModel_->kappa(alpha, Theta_, gs0_, rhoa_, da_, e_);

    // particle viscosity (Table 3.2, p.47)
    mua_ = viscosityModel_->mua(alpha, Theta_, gs0_, rhoa_, da_, e_);

    // dissipation (Eq. 3.24, p.50)
    volScalarField gammaCoeff =
        12.0*(1.0 - e_*e_)*sqr(alpha)*rhoa_*gs0_*(1.0/da_)
       *ThetaSqrt/sqrtPi;

    // Eq. 3.25, p. 50 Js = J1 - J2
    volScalarField J1 = 3.0*betaPrim;
    volScalarField J2 =
        0.25*sqr(betaPrim)*da_*sqr(Ur)
       /(alpha*rhoa_*sqrtPi*(ThetaSqrt + TsmallSqrt));

    // bulk viscosity  p. 45 (Lun et al. 1984).
    lambda_ = (4.0/3.0)*sqr(alpha_)*rhoa_*da_*gs0_*(1.0+e_)*ThetaSqrt/sqrtPi;


    // stress tensor, Definitions, Table 3.1, p. 43
    volTensorField tau = 2.0*mua_*D + (lambda_ - (2.0/3.0)*mua_)*tr(D)*I;

    if (!equilibrium_)
    {
        // construct the granular temperature equation (Eq. 3.20, p. 44)
        // NB. note that there are two typos in Eq. 3.20
        // no grad infront of Ps
        // wrong sign infront of laplacian
        fvScalarMatrix ThetaEqn
        (
            fvm::ddt(1.5*alpha*rhoa_, Theta_)
          + fvm::div(phi, Theta_, scheme)
         ==
            fvm::SuSp(-((PsCoeff*I) && dU), Theta_)
          + (tau && dU)
          + fvm::laplacian(kappa_, Theta_, "laplacian(kappa,Theta)")
          + fvm::Sp(-gammaCoeff, Theta_)
          + fvm::Sp(-J1, Theta_)
          + fvm::Sp(J2/(Theta_ + Tsmall), Theta_)
        );

        ThetaEqn.relax();
        ThetaEqn.solve();
    }
    else
    {
        // equilibrium => dissipation == production
        // Eq. 4.14, p.82
        volScalarField K1 = 2.0*(1.0 + e_)*rhoa_*gs0_;
        volScalarField K3 = 0.5*da_*rhoa_*
            (
                (sqrtPi/(3.0*(3.0-e_)))
               *(1.0 + 0.4*(1.0 + e_)*(3.0*e_ - 1.0)*alpha*gs0_)
              + 1.6*alpha*gs0_*(1.0 + e_)/sqrtPi
            );

        volScalarField K2 =
            4.0*da_*rhoa_*(1.0 + e_)*alpha*gs0_/(3.0*sqrtPi) - 2.0*K3/3.0;

        volScalarField K4 = 12.0*(1.0 - e_*e_)*rhoa_*gs0_/(da_*sqrtPi);

        volScalarField trD = tr(D);
        volTensorField D2 = D & D;
        volScalarField tr2D = trD*trD;
        volScalarField trD2 = tr(D2);

        volScalarField t1 = K1*alpha + rhoa_;
        volScalarField l1 = -t1*trD;
        volScalarField l2 = sqr(t1)*tr2D;
        volScalarField l3 = 4.0*K4*alpha*(2.0*K3*trD2 + K2*tr2D);

        Theta_ = sqr((l1 + sqrt(l2 + l3))/(2.0*(alpha + 1.0e-4)*K4));
    }

    Theta_.max(1.0e-15);
    Theta_.min(1.0e+3);

    volScalarField pf =
        frictionalStressModel_->frictionalPressure
    (
        alpha_,
        alphaMinFriction_,
        alphaMax_,
        Fr_,
        eta_,
        p_
    );

    PsCoeff += pf/(Theta_+Tsmall);

    PsCoeff.min(1.0e+10);
    PsCoeff.max(-1.0e+10);

    // update particle pressure
    pa_ = PsCoeff*Theta_;

    // frictional shear stress, Eq. 3.30, p. 52
    volScalarField muf = frictionalStressModel_->muf
    (
        alpha_,
        alphaMax_,
        pf,
        D,
        phi_
    );

    // add frictional stress for alpha > alphaMinFriction
    mua_ = viscosityModel_->mua(alpha, Theta_, gs0_, rhoa_, da_, e_) + muf;
    mua_.min(1.0e+2);
    mua_.max(0.0);

    lambda_ = (4.0/3.0)*sqr(alpha_)*rhoa_*da_*gs0_*(1.0 + e_)*ThetaSqrt/sqrtPi;

    Info<< "kinTheory: max(Theta) = " << max(Theta_).value() << endl;

    volScalarField ktn = mua_/rhoa_;

    Info<< "kinTheory: min(nua) = " << min(ktn).value()
        << ", max(nua) = " << max(ktn).value() << endl;

    Info<< "kinTheory: min(pa) = " << min(pa_).value()
        << ", max(pa) = " << max(pa_).value() << endl;
}


// ************************************************************************* //
