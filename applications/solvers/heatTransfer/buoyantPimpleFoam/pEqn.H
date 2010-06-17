{
    rho = thermo.rho();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    thermo.rho() -= psi*p_rgh;

    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rhorUAf("(rho*(1|A(U)))", fvc::interpolate(rho*rUA));

    U = rUA*UEqn.H();

    phi = fvc::interpolate(rho)*
    (
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

    surfaceScalarField buoyancyPhi = -rhorUAf*ghf*fvc::snGrad(rho)*mesh.magSf();
    phi += buoyancyPhi;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
          + fvc::div(phi)
          - fvm::laplacian(rhorUAf, p_rgh)
        );

        p_rghEqn.solve
        (
            mesh.solver
            (
                p_rgh.select
                (
                    (
                        finalIter
                     && corr == nCorr-1
                     && nonOrth == nNonOrthCorr
                    )
                )
            )
        );

        if (nonOrth == nNonOrthCorr)
        {
            // Calculate the conservative fluxes
            phi += p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U += rUA*fvc::reconstruct((buoyancyPhi + p_rghEqn.flux())/rhorUAf);
            U.correctBoundaryConditions();
        }
    }

    p = p_rgh + rho*gh;

    // Second part of thermodynamic density update
    thermo.rho() += psi*p_rgh;

    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"
}
