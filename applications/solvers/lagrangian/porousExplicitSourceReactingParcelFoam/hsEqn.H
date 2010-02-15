{
    tmp<volScalarField> pWork
    (
        new volScalarField
        (
            IOobject
            (
                "pWork",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );

    if (dpdt)
    {
        pWork() += fvc::ddt(p);
    }
    if (eWork)
    {
        pWork() = -p*fvc::div(phi/fvc::interpolate(rho));
    }
    if (hWork)
    {
        pWork() += fvc::div(phi/fvc::interpolate(rho)*fvc::interpolate(p));
    }

    {
        solve
        (
            fvm::ddt(rho, hs)
          + mvConvection->fvmDiv(phi, hs)
          - fvm::laplacian(turbulence->alphaEff(), hs)
         ==
            pWork()
          + parcels.Sh()
          + radiation->Shs(thermo)
          + energySource.Su()
          + chemistrySh
        );

        thermo.correct();

        radiation->correct();

        Info<< "T gas min/max   = " << min(T).value() << ", "
            << max(T).value() << endl;
    }
}
