fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
);

UEqn.relax();

solve
(
    UEqn
  ==
    fvc::reconstruct
    (
        (
            - ghf*fvc::snGrad(rho)
            - fvc::snGrad(p_rgh)
        )*mesh.magSf()
    ),
    mesh.solver(U.select(oCorr == nOuterCorr - 1))
);
