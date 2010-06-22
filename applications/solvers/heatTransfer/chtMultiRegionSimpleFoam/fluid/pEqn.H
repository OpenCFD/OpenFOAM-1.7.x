{
    rho = thermo.rho();
    rho = max(rho, rhoMin[i]);
    rho = min(rho, rhoMax[i]);
    rho.relax();

    volScalarField rUA = 1.0/UEqn().A();
    surfaceScalarField rhorUAf("(rho*(1|A(U)))", fvc::interpolate(rho*rUA));

    U = rUA*UEqn().H();
    UEqn.clear();

    phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());
    bool closedVolume = adjustPhi(phi, U, p_rgh);

    surfaceScalarField buoyancyPhi = rhorUAf*ghf*fvc::snGrad(rho)*mesh.magSf();
    phi -= buoyancyPhi;

    // Solve pressure
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rhorUAf, p_rgh) == fvc::div(phi)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        eqnResidual = p_rghEqn.solve().initialResidual();

        // Retain the residual from the first iteration
        if (nonOrth == 0)
        {
            maxResidual = max(eqnResidual, maxResidual);
        }

        if (nonOrth == nNonOrthCorr)
        {
            // Calculate the conservative fluxes
            phi -= p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U -= rUA*fvc::reconstruct((buoyancyPhi + p_rghEqn.flux())/rhorUAf);
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"

    p = p_rgh + rho*gh;

    // For closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (closedVolume)
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /fvc::domainIntegrate(psi);
        p_rgh = p - rho*gh;
    }

    rho = thermo.rho();
    rho = max(rho, rhoMin[i]);
    rho = min(rho, rhoMax[i]);
    rho.relax();

    Info<< "Min/max rho:" << min(rho).value() << ' '
        << max(rho).value() << endl;

    // Update thermal conductivity
    K = thermo.Cp()*turb.alphaEff();
}
