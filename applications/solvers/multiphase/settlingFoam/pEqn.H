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
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(rUAf, p_rgh) == fvc::ddt(rho) + fvc::div(phi)
    );

    p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));
    p_rghEqn.solve();

    if (nonOrth == nNonOrthCorr)
    {
        phi -= p_rghEqn.flux();
    }
}

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

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
U.correctBoundaryConditions();
