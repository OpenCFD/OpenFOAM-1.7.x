rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();

volScalarField rUA = 1.0/UEqn().A();
U = rUA*UEqn().H();
UEqn.clear();

bool closedVolume = false;

if (transonic)
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)*(fvc::interpolate(U) & mesh.Sf())
    );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::div(phid, p)
          - fvm::laplacian(rho*rUA, p)
        );

        // Relax the pressure equation to ensure diagonal-dominance
        pEqn.relax(mesh.relaxationFactor("pEqn"));

        pEqn.setReference(pRefCell, pRefValue);

        // retain the residual from the first iteration
        if (nonOrth == 0)
        {
            eqnResidual = pEqn.solve().initialResidual();
            maxResidual = max(eqnResidual, maxResidual);
        }
        else
        {
            pEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());
    closedVolume = adjustPhi(phi, U, p);

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rho*rUA, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        // Retain the residual from the first iteration
        if (nonOrth == 0)
        {
            eqnResidual = pEqn.solve().initialResidual();
            maxResidual = max(eqnResidual, maxResidual);
        }
        else
        {
            pEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pEqn.flux();
        }
    }
}


#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi*p))
        /fvc::domainIntegrate(psi);
}

rho = thermo.rho();
rho = max(rho, rhoMin);
rho = min(rho, rhoMax);
rho.relax();
Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;
