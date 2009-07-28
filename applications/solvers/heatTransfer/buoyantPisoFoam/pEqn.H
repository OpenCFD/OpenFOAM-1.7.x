{
    bool closedVolume = p.needReference();

    rho = thermo.rho();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    thermo.rho() -= psi*p;

    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rhorUAf("(rho*(1|A(U)))", fvc::interpolate(rho*rUA));

    U = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        )
    );

    phi = phiU + rhorUAf*fvc::interpolate(rho)*(g & mesh.Sf());

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phi)
          - fvm::laplacian(rhorUAf, p)
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
            phi += pEqn.flux();
        }
    }

    // Second part of thermodynamic density update
    thermo.rho() += psi*p;

    U += rUA*fvc::reconstruct((phi - phiU)/rhorUAf);
    U.correctBoundaryConditions();

    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        p +=
            (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
        thermo.rho() = psi*p;
        rho += (initialMass - fvc::domainIntegrate(rho))/totalVolume;
    }
}
