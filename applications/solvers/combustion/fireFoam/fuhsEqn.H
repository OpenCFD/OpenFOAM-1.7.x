{
    // Solve fuel equation
    // ~~~~~~~~~~~~~~~~~~~
    volScalarField fres = composition.fres(ft, stoicRatio.value());
    volScalarField wFuelNorm = combustion->wFuelNorm();

    {
        fvScalarMatrix fuEqn
        (
            fvm::ddt(rho, fu)
          + mvConvection->fvmDiv(phi, fu)
          - fvm::laplacian(turbulence->alphaEff(), fu)
         ==
          - fvm::Sp(wFuelNorm, fu)
          + wFuelNorm*fres
        );

        fuEqn.relax();
        fuEqn.solve();
    }

    volScalarField Rfu = wFuelNorm*(fu - fres);

    Info<< "max(fu) = " << max(fu).value() << endl;
    Info<< "min(fu) = " << min(fu).value() << endl;


    // Solve sensible enthalpy equation
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    radiation->correct();
    dQ = combustion->Sh(Rfu);

    {
        fvScalarMatrix hsEqn
        (
            fvm::ddt(rho, hs)
          + mvConvection->fvmDiv(phi,hs)
          - fvm::laplacian(turbulence->alphaEff(), hs)
         ==
            dpdt
          + fvc::div(phi/fvc::interpolate(rho)*fvc::interpolate(p))
          - p*fvc::div(phi/fvc::interpolate(rho))
          + dQ
          + radiation->Shs(thermo)
        );

        hsEqn.relax();
        hsEqn.solve();
    }

    thermo.correct();
    combustion->correct();
}
