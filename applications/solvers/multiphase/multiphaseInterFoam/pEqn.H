{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    U = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf()) + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

    adjustPhi(phiU, U, p_rgh);

    phi = phiU +
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*mesh.magSf()*rUAf;

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rUAf, p_rgh) == fvc::div(phi)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve
        (
            mesh.solver
            (
                p_rgh.select(corr == nCorr-1 && nonOrth == nNonOrthCorr)
            )
        );

        if (nonOrth == nNonOrthCorr)
        {
            phi -= p_rghEqn.flux();
        }
    }

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();

    #include "continuityErrs.H"

    p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }
}
