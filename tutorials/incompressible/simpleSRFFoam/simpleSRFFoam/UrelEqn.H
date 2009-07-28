    // Relative momentum predictor
    tmp<fvVectorMatrix> UrelEqn
    (
        fvm::div(phi, Urel)
      + turbulence->divDevReff(Urel)
      + SRF->Su()
    );

    UrelEqn().relax();

    eqnResidual = solve
    (
        UrelEqn() == -fvc::grad(p)
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);
