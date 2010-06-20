{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    tmp<fvScalarMatrix> p_rghEqnComp;

    if (transonic)
    {
        p_rghEqnComp =
        (
            fvm::ddt(p_rgh)
          + fvm::div(phi, p_rgh)
          - fvm::Sp(fvc::div(phi), p_rgh)
        );
    }
    else
    {
        p_rghEqnComp =
        (
            fvm::ddt(p_rgh)
          + fvc::div(phi, p_rgh)
          - fvc::Sp(fvc::div(phi), p_rgh)
        );
    }


    U = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

    phi = phiU +
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rUAf*mesh.magSf();

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqnIncomp
        (
            fvc::div(phi)
          - fvm::laplacian(rUAf, p_rgh)
        );

        solve
        (
            (
                max(alpha1, scalar(0))*(psi1/rho1)
              + max(alpha2, scalar(0))*(psi2/rho2)
            )
           *p_rghEqnComp()
          + p_rghEqnIncomp,
            mesh.solver
            (
                p_rgh.select
                (
                    oCorr == nOuterCorr-1
                 && corr == nCorr-1
                 && nonOrth == nNonOrthCorr
                )
            )
        );

        if (nonOrth == nNonOrthCorr)
        {
            dgdt =
                (pos(alpha2)*(psi2/rho2) - pos(alpha1)*(psi1/rho1))
               *(p_rghEqnComp & p_rgh);
            phi += p_rghEqnIncomp.flux();
        }
    }

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();

    p = max
    (
        (p_rgh + gh*(alpha1*rho10 + alpha2*rho20))
       /(1.0 - gh*(alpha1*psi1 + alpha2*psi2)),
        pMin
    );

    rho1 = rho10 + psi1*p;
    rho2 = rho20 + psi2*p;

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(p_rgh) " << min(p_rgh).value() << endl;
}
