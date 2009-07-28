{
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + mvConvection->fvmDiv(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        DpDt
     +  parcels.Sh()
     +  radiation->Sh(thermo)
    );

    hEqn.relax();

    hEqn.solve();

    thermo.correct();

    radiation->correct();
}
