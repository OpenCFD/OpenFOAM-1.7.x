{
    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
     ==
        DpDt
     +  coalParcels.Sh()
     +  limestoneParcels.Sh()
     +  enthalpySource.Su()
     +  radiation->Shs(thermo)
     +  chemistrySh
    );

    hsEqn.relax();

    hsEqn.solve();

    thermo.correct();

    radiation->correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
