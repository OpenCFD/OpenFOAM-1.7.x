    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    // Add resistance on the actuation disks
    actuationDisks.addSu(UEqn());

    UEqn().relax();

    eqnResidual = solve
    (
        UEqn() == -fvc::grad(p)
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);
