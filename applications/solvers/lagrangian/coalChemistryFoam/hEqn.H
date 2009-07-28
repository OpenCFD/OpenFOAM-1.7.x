{
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + mvConvection->fvmDiv(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        DpDt
     +  coalParcels.Sh()
     +  limestoneParcels.Sh()
     +  enthalpySource.Su()
     +  radiation->Sh(thermo)
    );

    hEqn.relax();

    hEqn.solve();

    thermo.correct();

    radiation->correct();
}
