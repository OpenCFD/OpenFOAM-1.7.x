{
    surfaceScalarField alphaf = fvc::interpolate(alpha);
    surfaceScalarField betaf = scalar(1) - alphaf;

    volScalarField rUaA = 1.0/UaEqn.A();
    volScalarField rUbA = 1.0/UbEqn.A();

    surfaceScalarField rUaAf = fvc::interpolate(rUaA);
    surfaceScalarField rUbAf = fvc::interpolate(rUbA);

    Ua = rUaA*UaEqn.H();
    Ub = rUbA*UbEqn.H();

    surfaceScalarField phiDraga =
        fvc::interpolate(beta/rhoa*dragCoef*rUaA)*phib + rUaAf*(g & mesh.Sf());
    surfaceScalarField phiDragb =
        fvc::interpolate(alpha/rhob*dragCoef*rUbA)*phia + rUbAf*(g & mesh.Sf());

    forAll(p.boundaryField(), patchi)
    {
        if (isA<zeroGradientFvPatchScalarField>(p.boundaryField()[patchi]))
        {
            phiDraga.boundaryField()[patchi] = 0.0;
            phiDragb.boundaryField()[patchi] = 0.0;
        }
    }

    phia = (fvc::interpolate(Ua) & mesh.Sf()) + fvc::ddtPhiCorr(rUaA, Ua, phia)
        + phiDraga;
    phib = (fvc::interpolate(Ub) & mesh.Sf()) + fvc::ddtPhiCorr(rUbA, Ub, phib)
        + phiDragb;

    phi = alphaf*phia + betaf*phib;

    surfaceScalarField Dp("(rho*(1|A(U)))", alphaf*rUaAf/rhoa + betaf*rUbAf/rhob);

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(Dp, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            surfaceScalarField SfGradp = pEqn.flux()/Dp;

            phia -= rUaAf*SfGradp/rhoa;
            phib -= rUbAf*SfGradp/rhob;
            phi = alphaf*phia + betaf*phib;

            p.relax();
            SfGradp = pEqn.flux()/Dp;

            Ua += (fvc::reconstruct(phiDraga - rUaAf*SfGradp/rhoa));
            //Ua += rUaA*(fvc::reconstruct(phiDraga/rUaAf - SfGradp/rhoa));
            Ua.correctBoundaryConditions();

            Ub += (fvc::reconstruct(phiDragb - rUbAf*SfGradp/rhob));
            //Ub += rUbA*(fvc::reconstruct(phiDragb/rUbAf - SfGradp/rhob));
            Ub.correctBoundaryConditions();

            U = alpha*Ua + beta*Ub;
        }
    }
}

#include "continuityErrs.H"
