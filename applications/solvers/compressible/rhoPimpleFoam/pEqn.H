rho = thermo.rho();

volScalarField rUA = 1.0/UEqn().A();
U = rUA*UEqn().H();

if (nCorr <= 1)
{
    UEqn.clear();
}

if (transonic)
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rUA, rho, U, phi)
        )
    );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*rUA, p)
        );

        if
        (
            oCorr == nOuterCorr-1
            && corr == nCorr-1
            && nonOrth == nNonOrthCorr
        )
        {
            pEqn.solve(mesh.solver("pFinal"));
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
    phi =
        fvc::interpolate(rho)*
        (
            (fvc::interpolate(U) & mesh.Sf())
        //+ fvc::ddtPhiCorr(rUA, rho, U, phi)
        );

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*rUA, p)
        );

        if
        (
            oCorr == nOuterCorr-1
         && corr == nCorr-1
         && nonOrth == nNonOrthCorr
        )
        {
            pEqn.solve(mesh.solver("pFinal"));
        }
        else
        {
            pEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

//if (oCorr != nOuterCorr-1)
{
    // Explicitly relax pressure for momentum corrector
    p.relax();

    rho = thermo.rho();
    rho.relax();
    Info<< "rho max/min : " << max(rho).value()
        << " " << min(rho).value() << endl;
}

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

bound(p, pMin);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
/*
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi*p))
        /fvc::domainIntegrate(psi);
}
*/
