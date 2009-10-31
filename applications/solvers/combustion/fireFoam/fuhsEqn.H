{
    // Solve fuel equation
    // ~~~~~~~~~~~~~~~~~~~
    fvScalarMatrix R = combustion->R(fu);

    {
        fvScalarMatrix fuEqn
        (
            fvm::ddt(rho, fu)
          + mvConvection->fvmDiv(phi, fu)
          - fvm::laplacian(turbulence->alphaEff(), fu)
         ==
            R
        );

        fuEqn.relax();
        fuEqn.solve();
    }

    Info<< "max(fu) = " << max(fu).value() << endl;
    Info<< "min(fu) = " << min(fu).value() << endl;


    // Solve sensible enthalpy equation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    radiation->correct();
    dQ = combustion->dQ(R);

    {
        fvScalarMatrix hsEqn
        (
            fvm::ddt(rho, hs)
          + mvConvection->fvmDiv(phi,hs)
          - fvm::laplacian(turbulence->alphaEff(), hs)
         ==
            DpDt
          + dQ
          + radiation->Shs(thermo)
        );

        hsEqn.relax();
        hsEqn.solve();
    }

    thermo.correct();
    combustion->correct();
}
