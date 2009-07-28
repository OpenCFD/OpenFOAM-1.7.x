tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,ft_b_h_hu)")
    )
);

if (composition.contains("ft"))
{
    volScalarField& ft = composition.Y("ft");

    solve
    (
        fvm::ddt(rho, ft)
      + mvConvection->fvmDiv(phi, ft)
      - fvm::laplacian(turbulence->alphaEff(), ft)
    );
}
