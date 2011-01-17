volScalarField rAUrel = 1.0/UrelEqn().A();
Urel = rAUrel*UrelEqn().H();

if (nCorr <= 1)
{
    UrelEqn.clear();
}

phi = (fvc::interpolate(Urel) & mesh.Sf())
    + fvc::ddtPhiCorr(rAUrel, Urel, phi);

adjustPhi(phi, Urel, p);

// Non-orthogonal pressure corrector loop
for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAUrel, p) == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);

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
        phi -= pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector except for last corrector
if (oCorr != nOuterCorr-1)
{
    p.relax();
}

// Momentum corrector
Urel -= rAUrel*fvc::grad(p);
Urel.correctBoundaryConditions();
