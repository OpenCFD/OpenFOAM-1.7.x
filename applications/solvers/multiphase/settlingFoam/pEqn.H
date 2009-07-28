volScalarField rUA = 1.0/UEqn.A();

surfaceScalarField rUAf
(
    "(rho*(1|A(U)))",
    fvc::interpolate(rho)*fvc::interpolate(rUA)
);

U = rUA*UEqn.H();
phi =
    fvc::interpolate(rho)
   *(
       (fvc::interpolate(U) & mesh.Sf())
     + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

surfaceScalarField phiU("phiU", phi);
phi += fvc::interpolate(rho)*(g & mesh.Sf())*rUAf;

for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rUAf, p) == fvc::ddt(rho) + fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve();

    if (nonOrth == nNonOrthCorr)
    {
        phi -= pEqn.flux();
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
U.correctBoundaryConditions();
