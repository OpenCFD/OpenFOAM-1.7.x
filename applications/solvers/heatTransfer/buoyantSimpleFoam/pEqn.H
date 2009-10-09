{
    rho = thermo.rho();

    volScalarField rUA = 1.0/UEqn().A();
    surfaceScalarField rhorUAf("(rho*(1|A(U)))", fvc::interpolate(rho*rUA));

    U = rUA*UEqn().H();
    UEqn.clear();

    phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());
    bool closedVolume = adjustPhi(phi, U, p);

    surfaceScalarField buoyancyPhi =
        rhorUAf*fvc::interpolate(rho)*(g & mesh.Sf());
    phi += buoyancyPhi;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rhorUAf, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        // retain the residual from the first iteration
        if (nonOrth == 0)
        {
            eqnResidual = pEqn.solve().initialResidual();
            maxResidual = max(eqnResidual, maxResidual);
        }
        else
        {
            pEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            // For closed-volume cases adjust the pressure and density levels
            // to obey overall mass continuity
            if (closedVolume)
            {
                p += (initialMass - fvc::domainIntegrate(psi*p))
                    /fvc::domainIntegrate(psi);
            }

            // Calculate the conservative fluxes
            phi -= pEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U += rUA*fvc::reconstruct((buoyancyPhi - pEqn.flux())/rhorUAf);
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"

    rho = thermo.rho();
    rho.relax();
    Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value()
        << endl;
}
