{
    solve
    (
        betav*fvm::ddt(rho, h)
      + mvConvection->fvmDiv(phi, h)
      - fvm::laplacian(Db, h)
     ==
        betav*DpDt
    );

    thermo.correct();
}
