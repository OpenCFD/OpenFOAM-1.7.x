bool closedVolume = false;

rho = thermo.rho();

volScalarField rUA = 1.0/UEqn.A();
surfaceScalarField rhorUAf("(rho*(1|A(U)))", fvc::interpolate(rho*rUA));
U = rUA*UEqn.H();

surfaceScalarField phiU
(
    fvc::interpolate(rho)
   *(
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho, U, phi)
    )
);

phi = phiU - rhorUAf*ghf*fvc::snGrad(rho)*mesh.magSf();

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    surfaceScalarField rhorUAf = fvc::interpolate(rho*rUA);

    fvScalarMatrix p_rghEqn
    (
        fvm::ddt(psi, p_rgh)
      + fvc::div(phi)
      - fvm::laplacian(rhorUAf, p_rgh)
    );

    closedVolume = p_rgh.needReference();

    p_rghEqn.solve
    (
        mesh.solver(p_rgh.select(corr == nCorr-1 && nonOrth == nNonOrthCorr))
    );

    if (nonOrth == nNonOrthCorr)
    {
        phi += p_rghEqn.flux();
    }
}

p = p_rgh + rho*gh;

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U += rUA*fvc::reconstruct((phi - phiU)/rhorUAf);
U.correctBoundaryConditions();

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p +=
        (initialMass - fvc::domainIntegrate(thermo.psi()*p))
       /fvc::domainIntegrate(thermo.psi());
    p_rgh == p - rho*gh;
    rho = thermo.rho();
}

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
