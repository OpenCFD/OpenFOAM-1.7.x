{
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + fvm::div(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        DpDt
      + thermoCloud1.Sh()
    );

    hEqn.relax();

    hEqn.solve();

    thermo.correct();
}
