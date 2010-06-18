    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
    );

    UEqn().relax();

    if (momentumPredictor)
    {
        eqnResidual = solve
        (
            UEqn()
         ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        ).initialResidual();

        maxResidual = max(eqnResidual, maxResidual);
    }
