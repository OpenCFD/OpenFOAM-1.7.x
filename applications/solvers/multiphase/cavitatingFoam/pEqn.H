{
    if (nOuterCorr == 1)
    {
        p =
        (
            rho
          - (1.0 - gamma)*rhol0
          - ((gamma*psiv + (1.0 - gamma)*psil) - psi)*pSat
        )/psi;
    }

    surfaceScalarField rhof = fvc::interpolate(rho, "rhof");

    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf("rUAf", rhof*fvc::interpolate(rUA));
    volVectorField HbyA = rUA*UEqn.H();

    phiv = (fvc::interpolate(HbyA) & mesh.Sf())
         + fvc::ddtPhiCorr(rUA, rho, U, phiv);

    p.boundaryField().updateCoeffs();

    surfaceScalarField phiGradp = rUAf*mesh.magSf()*fvc::snGrad(p);

    phiv -= phiGradp/rhof;

    #include "resetPhivPatches.H"

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          - (rhol0 + (psil - psiv)*pSat)*fvc::ddt(gamma) - pSat*fvc::ddt(psi)
          + fvc::div(phiv, rho)
          + fvc::div(phiGradp)
          - fvm::laplacian(rUAf, p)
        );

        if (corr == nCorr-1 && nonOrth == nNonOrthCorr)
        {
            pEqn.solve(mesh.solver(p.name() + "Final"));
        }
        else
        {
            pEqn.solve(mesh.solver(p.name()));
        }

        if (nonOrth == nNonOrthCorr)
        {
            phiv += (phiGradp + pEqn.flux())/rhof;
        }
    }

    Info<< "Predicted p max-min : " << max(p).value()
        << " " << min(p).value() << endl;

    rho == max
    (
        psi*p
      + (1.0 - gamma)*rhol0
      + ((gamma*psiv + (1.0 - gamma)*psil) - psi)*pSat,
        rhoMin
    );

    #include "gammaPsi.H"

    p =
    (
        rho
      - (1.0 - gamma)*rhol0
      - ((gamma*psiv + (1.0 - gamma)*psil) - psi)*pSat
    )/psi;

    p.correctBoundaryConditions();

    Info<< "Phase-change corrected p max-min : " << max(p).value()
        << " " << min(p).value() << endl;

    // Correct velocity

    U = HbyA - rUA*fvc::grad(p);

    // Remove the swirl component of velocity for "wedge" cases
    if (piso.found("removeSwirl"))
    {
        label swirlCmpt(readLabel(piso.lookup("removeSwirl")));

        Info<< "Removing swirl component-" << swirlCmpt << " of U" << endl;
        U.field().replace(swirlCmpt, 0.0);
    }

    U.correctBoundaryConditions();

    Info<< "max(U) " << max(mag(U)).value() << endl;
}
