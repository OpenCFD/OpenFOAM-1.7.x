{
    volScalarField kappaEff
    (
        "kappaEff",
        turbulence->nu()/Pr + turbulence->nut()/Prt
    );

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
    );

    TEqn.relax();
    TEqn.solve(mesh.solver(T.select(finalIter)));

    rhok = 1.0 - beta*(T - TRef);
}
