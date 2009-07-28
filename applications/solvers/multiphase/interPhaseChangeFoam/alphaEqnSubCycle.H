surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0)
);

{
    label nAlphaCorr
    (
        readLabel(piso.lookup("nAlphaCorr"))
    );

    label nAlphaSubCycles
    (
        readLabel(piso.lookup("nAlphaSubCycles"))
    );

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(interface.cAlpha()*phic, max(phic));

    volScalarField divU = fvc::div(phi);

    dimensionedScalar totalDeltaT = runTime.deltaT();

    if (nAlphaSubCycles > 1)
    {
        for
        (
            subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
            !(++alphaSubCycle).end();
        )
        {
#           include "alphaEqn.H"
        }
    }
    else
    {
#       include "alphaEqn.H"
    }

    if (nOuterCorr == 1)
    {
        interface.correct();
    }

    rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;
}
