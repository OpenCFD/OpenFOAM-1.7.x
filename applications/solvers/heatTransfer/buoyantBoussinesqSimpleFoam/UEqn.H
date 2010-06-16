    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
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
                  - ghf*fvc::snGrad(rhok)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        ).initialResidual();

        maxResidual = max(eqnResidual, maxResidual);
    }
