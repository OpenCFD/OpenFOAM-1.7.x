tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,ft_b_h)")
    )
);

{
    fvScalarMatrix ftEqn
    (
        fvm::ddt(rho, ft)
      + mvConvection->fvmDiv(phi, ft)
      - fvm::laplacian(turbulence->alphaEff(), ft)
    );

    ftEqn.relax();
    ftEqn.solve();
}

Info<< "max(ft) = " << max(ft).value() << endl;
Info<< "min(ft) = " << min(ft).value() << endl;
