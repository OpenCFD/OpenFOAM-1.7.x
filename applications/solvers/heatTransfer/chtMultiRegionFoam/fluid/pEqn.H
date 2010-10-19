{
    bool closedVolume = p_rgh.needReference();
    dimensionedScalar compressibility = fvc::domainIntegrate(psi);
    bool compressible = (compressibility.value() > SMALL);

    rho = thermo.rho();

    volScalarField rUA = 1.0/UEqn().A();
    surfaceScalarField rhorUAf("(rho*(1|A(U)))", fvc::interpolate(rho*rUA));

    U = rUA*UEqn().H();

    surfaceScalarField phiU
    (
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        )
    );

    phi = phiU - rhorUAf*ghf*fvc::snGrad(rho)*mesh.magSf();

    {
        fvScalarMatrix p_rghDDtEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
          + fvc::div(phi)
        );

        // Thermodynamic density needs to be updated by psi*d(p) after the
        // pressure solution - done in 2 parts. Part 1:
        thermo.rho() -= psi*p_rgh;

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix p_rghEqn
            (
                p_rghDDtEqn
              - fvm::laplacian(rhorUAf, p_rgh)
            );

            p_rghEqn.solve
            (
                mesh.solver
                (
                    p_rgh.select
                    (
                        (
                           oCorr == nOuterCorr-1
                        && corr == nCorr-1
                        && nonOrth == nNonOrthCorr
                        )
                    )
                )
            );

            if (nonOrth == nNonOrthCorr)
            {
                phi += p_rghEqn.flux();
            }
        }

        // Second part of thermodynamic density update
        thermo.rho() += psi*p_rgh;
    }

    // Correct velocity field
    U += rUA*fvc::reconstruct((phi - phiU)/rhorUAf);
    U.correctBoundaryConditions();

    p = p_rgh + rho*gh;

    // Update pressure substantive derivative
    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

    // Solve continuity
    #include "rhoEqn.H"

    // Update continuity errors
    #include "compressibleContinuityErrors.H"

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume && compressible)
    {
        p += (initialMass - fvc::domainIntegrate(thermo.rho()))
            /compressibility;
        rho = thermo.rho();
        p_rgh = p - rho*gh;
    }

    // Update thermal conductivity
    K = thermoFluid[i].Cp()*turb.alphaEff();
}
