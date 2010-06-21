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

    if (finalIter)
    {
        UEqn.relax(1);
    }
    else
    {
        UEqn.relax();
    }

    if (momentumPredictor)
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            ),
            mesh.solver(U.select(finalIter))
        );
    }
