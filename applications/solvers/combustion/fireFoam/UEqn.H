fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
);

UEqn.relax();

if (oCorr == nOuterCorr - 1)
{
    solve
    (
        UEqn
      ==
        fvc::reconstruct
        (
            fvc::interpolate(rho)*(g & mesh.Sf())
          - fvc::snGrad(p)*mesh.magSf()
        ),
        mesh.solver("UFinal")
    );
}
else
{
    solve
    (
        UEqn
      ==
        fvc::reconstruct
        (
            fvc::interpolate(rho)*(g & mesh.Sf())
          - fvc::snGrad(p)*mesh.magSf()
        )
    );
}
