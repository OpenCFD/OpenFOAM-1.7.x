    fvVectorMatrix UEqn
    (
//        fvm::ddt(rho, U)
        pZones.ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho.dimensionedInternalField()*g
      + parcels.SU()
      + momentumSource.Su()
    );

    pZones.addResistance(UEqn);

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }

