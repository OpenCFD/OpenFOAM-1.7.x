{
    volScalarField rAU = 1.0/UEqn.A();
    surfaceScalarField rAUf = fvc::interpolate(rAU);

    U = rAU*UEqn.H();
    surfaceScalarField phiU("phiU", (fvc::interpolate(U) & mesh.Sf()));

    if (p.needReference())
    {
        fvc::makeRelative(phiU, U);
        adjustPhi(phiU, U, p);
        fvc::makeAbsolute(phiU, U);
    }

    phi = phiU +
    (
        fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)*mesh.magSf()
      + fvc::interpolate(rho)*(g & mesh.Sf())
    )*rAUf;

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUf, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

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
            phi -= pEqn.flux();
        }
    }

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();

    #include "continuityErrs.H"

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
}
