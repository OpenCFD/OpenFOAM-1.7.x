{
    p.boundaryField().updateCoeffs();
    volScalarField AUrel = UrelEqn().A();
    Urel = UrelEqn().H()/AUrel;
    UrelEqn.clear();
    phi = fvc::interpolate(Urel) & mesh.Sf();
    adjustPhi(phi, Urel, p);

    // Non-orthogonal pressure corrector loop
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(1.0/AUrel, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pEqn.flux();
       }
    }

#   include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    Urel -= fvc::grad(p)/AUrel;
    Urel.correctBoundaryConditions();
}
