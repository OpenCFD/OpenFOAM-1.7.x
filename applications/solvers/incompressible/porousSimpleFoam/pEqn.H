if (pressureImplicitPorosity)
{
    U = trTU()&UEqn().H();
}
else
{
    U = trAU()*UEqn().H();
}

UEqn.clear();
phi = fvc::interpolate(U) & mesh.Sf();
adjustPhi(phi, U, p);

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    tmp<fvScalarMatrix> tpEqn;

    if (pressureImplicitPorosity)
    {
        tpEqn = (fvm::laplacian(trTU(), p) == fvc::div(phi));
    }
    else
    {
        tpEqn = (fvm::laplacian(trAU(), p) == fvc::div(phi));
    }

    tpEqn().setReference(pRefCell, pRefValue);
    // retain the residual from the first iteration
    if (nonOrth == 0)
    {
        eqnResidual = tpEqn().solve().initialResidual();
        maxResidual = max(eqnResidual, maxResidual);
    }
    else
    {
        tpEqn().solve();
    }

    if (nonOrth == nNonOrthCorr)
    {
        phi -= tpEqn().flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

if (pressureImplicitPorosity)
{
    U -= trTU()&fvc::grad(p);
}
else
{
    U -= trAU()*fvc::grad(p);
}

U.correctBoundaryConditions();
