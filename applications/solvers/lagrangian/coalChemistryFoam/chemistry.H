{
    Info<< "Solving chemistry" << endl;

    chemistry.solve
    (
        runTime.value() - runTime.deltaT().value(),
        runTime.deltaT().value()
    );

    // turbulent time scale
    if (turbulentReaction)
    {
        DimensionedField<scalar, volMesh> tk =
            Cmix*sqrt(turbulence->muEff()/rho/turbulence->epsilon());
        DimensionedField<scalar, volMesh> tc =
            chemistry.tc()().dimensionedInternalField();

        // Chalmers PaSR model
        kappa = (runTime.deltaT() + tc)/(runTime.deltaT() + tc + tk);
    }
    else
    {
        kappa = 1.0;
    }

    chemistrySh = kappa*chemistry.Sh()();
}
