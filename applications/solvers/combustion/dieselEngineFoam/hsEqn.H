{
    solve
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
     ==
       DpDt
     + dieselSpray.heatTransferSource()().dimensionedInternalField()
     + chemistrySh
    );

    thermo.correct();
}
