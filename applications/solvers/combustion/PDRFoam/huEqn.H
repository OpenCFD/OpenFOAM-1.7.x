if (ign.ignited())
{
    solve
    (
        betav*fvm::ddt(rho, hu)
      + mvConvection->fvmDiv(phi, hu)
      - fvm::laplacian(Db, hu)

    // These terms cannot be used in partially-premixed combustion due to
    // the resultant inconsistency between ft and hu transport.
    // A possible solution would be to solve for ftu as well as ft.
    //- fvm::div(muEff*fvc::grad(b)/(b + 0.001), hu)
    //+ fvm::Sp(fvc::div(muEff*fvc::grad(b)/(b + 0.001)), hu)

     ==
        betav*DpDt*rho/thermo.rhou()
    );
}
