{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    U = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

    adjustPhi(phiU, U, p);

    phi = phiU +
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)*mesh.magSf()
          + fvc::interpolate(rho)*(g & mesh.Sf())
        )*rUAf;

    Pair<tmp<volScalarField> > vDotP = twoPhaseProperties->vDotP();
    const volScalarField& vDotcP = vDotP[0]();
    const volScalarField& vDotvP = vDotP[1]();

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvc::div(phi) - fvm::laplacian(rUAf, p)
          - (vDotvP - vDotcP)*pSat + fvm::Sp(vDotvP - vDotcP, p)
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
            phi += pEqn.flux();
        }
    }

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();
}
