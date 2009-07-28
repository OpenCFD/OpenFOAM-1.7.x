{
    fvScalarMatrix alpha1Eqn
    (
        fvm::ddt(alpha1)
      + fvm::div(phi, alpha1)
      //- fvm::Sp(fvc::div(phi), alpha1)
      - fvm::laplacian
        (
            Dab + alphatab*turbulence->nut(), alpha1,
            "laplacian(Dab,alpha1)"
        )
    );

    alpha1Eqn.solve();

    rhoPhi = alpha1Eqn.flux()*(rho1 - rho2) + phi*rho2;
    rho = alpha1*rho1 + (scalar(1) - alpha1)*rho2;

    Info<< "Phase 1 volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
