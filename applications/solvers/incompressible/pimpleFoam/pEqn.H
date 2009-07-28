U = rUA*UEqn().H();

if (nCorr <= 1)
{
    UEqn.clear();
}

phi = (fvc::interpolate(U) & mesh.Sf())
    + fvc::ddtPhiCorr(rUA, U, phi);

adjustPhi(phi, U, p);

// Non-orthogonal pressure corrector loop
for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rUA, p) == fvc::div(phi)
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

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();
