    // Solve the Momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + fvc::div
        (
            (Alpha/(scalar(1.001) - Alpha))*(sqr(rhoc)/rho)*Vdj*Vdj,
            "div(phiVdj,Vdj)"
        )
      - fvm::laplacian(mu, U, "laplacian(muEff,U)")
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
              - fvc::snGrad(p)*mesh.magSf()
            )
        );
    }
