{
    volScalarField rUA("rUA", 1.0/UEqn().A());
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

    U = rUA*UEqn().H();
    UEqn.clear();

    phi = fvc::interpolate(U) & mesh.Sf();
    adjustPhi(phi, U, p_rgh);

    surfaceScalarField buoyancyPhi = rUAf*ghf*fvc::snGrad(rhok)*mesh.magSf();
    phi -= buoyancyPhi;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rUAf, p_rgh) == fvc::div(phi)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        // retain the residual from the first iteration
        if (nonOrth == 0)
        {
            eqnResidual = p_rghEqn.solve().initialResidual();
            maxResidual = max(eqnResidual, maxResidual);
        }
        else
        {
            p_rghEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            // Calculate the conservative fluxes
            phi -= p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U -= rUA*fvc::reconstruct((buoyancyPhi + p_rghEqn.flux())/rUAf);
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"

    p = p_rgh + rhok*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rhok*gh;
    }
}
