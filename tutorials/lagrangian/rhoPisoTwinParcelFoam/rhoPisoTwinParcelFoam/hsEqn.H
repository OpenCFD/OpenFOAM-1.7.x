{
    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + fvm::div(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
     ==
        DpDt
      + thermoCloud1.Sh()
    );

    hsEqn.relax();

    hsEqn.solve();

    thermo.correct();
}
