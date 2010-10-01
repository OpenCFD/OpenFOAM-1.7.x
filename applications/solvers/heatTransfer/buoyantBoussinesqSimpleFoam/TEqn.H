{
    kappat = turbulence->nut()/Prt;
    kappat.correctBoundaryConditions();

    volScalarField kappaEff("kappaEff", turbulence->nu()/Pr + kappat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
    );

    TEqn.relax();

    eqnResidual = TEqn.solve().initialResidual();
    maxResidual = max(eqnResidual, maxResidual);

    rhok = 1.0 - beta*(T - TRef);
}
