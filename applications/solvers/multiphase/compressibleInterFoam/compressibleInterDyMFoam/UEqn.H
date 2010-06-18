    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U)
      - (fvc::grad(U) & fvc::grad(muEff))
    //- fvc::div(muf*(mesh.Sf() & fvc::interpolate(fvc::grad(U)().T())))
    );

    UEqn.relax();

    if (momentumPredictor)
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                fvc::interpolate(rho)*(g & mesh.Sf())
              + (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - fvc::snGrad(p)
                ) * mesh.magSf()
            )
        );
    }
