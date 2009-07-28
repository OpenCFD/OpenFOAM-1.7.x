{
    volScalarField kappaEff
    (
        "kappaEff",
        turbulence->nu()/Pr + turbulence->nut()/Prt
    );

    fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      - fvm::Sp(fvc::div(phi), T)
      - fvm::laplacian(kappaEff, T)
    );

    TEqn.relax();

    eqnResidual = TEqn.solve().initialResidual();
    maxResidual = max(eqnResidual, maxResidual);

    rhok = 1.0 - beta*(T - TRef);
}
