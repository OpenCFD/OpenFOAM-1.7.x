{
    solve
    (
        fvm::ddt(rho, e)
      + fvm::div(phi, e)
      - fvm::laplacian(turbulence->alphaEff(), e)
      ==
      - p*fvc::div(phi/fvc::interpolate(rho))
    );

    thermo.correct();
}
