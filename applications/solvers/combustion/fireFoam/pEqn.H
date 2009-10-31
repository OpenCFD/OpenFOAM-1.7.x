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

phi = phiU + rhorUAf*fvc::interpolate(rho)*(g & mesh.Sf());

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
       surfaceScalarField rhorUAf = fvc::interpolate(rho*rUA);

        fvScalarMatrix pEqn
        (
            fvm::ddt(psi,p)
          + fvc::div(phi)
          - fvm::laplacian(rhorUAf, p)
        );

        closedVolume = p.needReference();

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

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

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
    rho = thermo.rho();
}
