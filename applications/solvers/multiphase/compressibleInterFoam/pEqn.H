{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    tmp<fvScalarMatrix> pEqnComp;

    if (transonic)
    {
        pEqnComp =
            (fvm::ddt(p) + fvm::div(phi, p) - fvm::Sp(fvc::div(phi), p));
    }
    else
    {
        pEqnComp =
            (fvm::ddt(p) + fvc::div(phi, p) - fvc::Sp(fvc::div(phi), p));
    }


    U = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf()) + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

    phi = phiU +
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)*mesh.magSf()
          + fvc::interpolate(rho)*(g & mesh.Sf())
        )*rUAf;

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phi)
          - fvm::laplacian(rUAf, p)
        );

        solve
        (
            (
                max(alpha1, scalar(0))*(psi1/rho1)
              + max(alpha2, scalar(0))*(psi2/rho2)
            )
           *pEqnComp()
          + pEqnIncomp
        );

        if (nonOrth == nNonOrthCorr)
        {
            dgdt =
                (pos(alpha2)*(psi2/rho2) - pos(alpha1)*(psi1/rho1))
               *(pEqnComp & p);
            phi += pEqnIncomp.flux();
        }
    }

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();

    p.max(pMin);

    rho1 = rho10 + psi1*p;
    rho2 = rho20 + psi2*p;

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p) " << min(p).value() << endl;
}
