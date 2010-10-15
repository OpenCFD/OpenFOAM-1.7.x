{
    rho = thermo.rho();

    volScalarField rUA = 1.0/UEqn.A();
    U = rUA*UEqn.H();

    if (transonic)
    {
        surfaceScalarField phiv =
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi);

        phi = fvc::interpolate(rho)*phiv;

        surfaceScalarField phid
        (
            "phid",
            fvc::interpolate(thermo.psi())*phiv
        );

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + fvc::div(phi)
          + correction(fvm::ddt(psi, p) + fvm::div(phid, p))
        );

        // Thermodynamic density needs to be updated by psi*d(p) after the
        // pressure solution - done in 2 parts. Part 1:
        thermo.rho() -= psi*p;

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                pDDtEqn - fvm::laplacian(rho*rUA, p)
            );

            if (ocorr == nOuterCorr && corr == nCorr && nonOrth == nNonOrthCorr)
            {
                pEqn.solve(mesh.solver(p.name() + "Final"));
            }
            else
            {
                pEqn.solve();
            }

            if (nonOrth == nNonOrthCorr)
            {
                phi += pEqn.flux();
            }
        }

        // Second part of thermodynamic density update
        thermo.rho() += psi*p;
    }
    else
    {
        phi =
            fvc::interpolate(rho)
           *(
                (fvc::interpolate(U) & mesh.Sf())
              + fvc::ddtPhiCorr(rUA, rho, U, phi)
            );

        fvScalarMatrix pDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phi)
        );

        // Thermodynamic density needs to be updated by psi*d(p) after the
        // pressure solution - done in 2 parts. Part 1:
        thermo.rho() -= psi*p;

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                pDDtEqn - fvm::laplacian(rho*rUA, p)
            );

            if (ocorr == nOuterCorr && corr == nCorr && nonOrth == nNonOrthCorr)
            {
                pEqn.solve(mesh.solver(p.name() + "Final"));
            }
            else
            {
                pEqn.solve();
            }

            if (nonOrth == nNonOrthCorr)
            {
                phi += pEqn.flux();
            }
        }

        // Second part of thermodynamic density update
        thermo.rho() += psi*p;
    }

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"

    U -= rUA*fvc::grad(p);
    U.correctBoundaryConditions();

    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
}
