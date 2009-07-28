{
    bool closedVolume = p.needReference();

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

    phi = phiU + fvc::interpolate(rho)*(g & mesh.Sf())*rhorUAf;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rhorUAf, p)
        );

        if
        (
            oCorr == nOuterCorr-1
         && corr == nCorr-1
         && nonOrth == nNonOrthCorr
        )
        {
            pEqn.solve(mesh.solver(p.name() + "Final"));
        }
        else
        {
            pEqn.solve(mesh.solver(p.name()));
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
    }

    // Correct velocity field
    U += rUA*fvc::reconstruct((phi - phiU)/rhorUAf);
    U.correctBoundaryConditions();

    // Update pressure substantive derivative
    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

    // Solve continuity
    #include "rhoEqn.H"

    // Update continuity errors
    #include "compressibleContinuityErrors.H"

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (massIni - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
        rho = thermo.rho();
    }

    // Update thermal conductivity
    K = thermoFluid[i].Cp()*turb.alphaEff();
}
