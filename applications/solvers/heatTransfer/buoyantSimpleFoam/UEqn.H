    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
    );

    UEqn().relax();

    eqnResidual = solve
    (
        UEqn()
     ==
        fvc::reconstruct
        (
            fvc::interpolate(rho)*(g & mesh.Sf())
          - fvc::snGrad(p)*mesh.magSf()
        )
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);
