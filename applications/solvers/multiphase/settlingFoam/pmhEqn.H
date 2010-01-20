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
phi -= ghf*fvc::snGrad(rho)*rUAf*mesh.magSf();

for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pmhEqn
    (
        fvm::laplacian(rUAf, pmh) == fvc::ddt(rho) + fvc::div(phi)
    );

    pmhEqn.setReference(pmhRefCell, pmhRefValue);
    pmhEqn.solve();

    if (nonOrth == nNonOrthCorr)
    {
        phi -= pmhEqn.flux();
    }
}

p == pmh + rho*gh;

if (pmh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pmhRefCell)
    );
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
U.correctBoundaryConditions();
