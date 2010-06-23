    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    UEqn().relax();

    // Add resistance on the actuation disks
    actuationDisks.addSu(UEqn());

    eqnResidual = solve
    (
        UEqn() == -fvc::grad(p)
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);