    // Relative momentum predictor
    tmp<fvVectorMatrix> UrelEqn
    (
        fvm::ddt(Urel)
      + fvm::div(phi, Urel)
      + turbulence->divDevReff(Urel)
      + SRF->Su()
    );

    UrelEqn().relax();

    solve(UrelEqn() == -fvc::grad(p));
