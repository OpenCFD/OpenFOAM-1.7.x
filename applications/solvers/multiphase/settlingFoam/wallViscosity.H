{
    scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    scalar kappa_ = kappa.value();
    scalar E_ = E.value();

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isType<wallFvPatch>(curPatch))
        {
            const scalarField& rhow = rho.boundaryField()[patchi];

            const scalarField muw = mul.boundaryField()[patchi];
            scalarField& mutw = mut.boundaryField()[patchi];

            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                scalar yPlus =
                    Cmu25*y[patchi][facei]*::sqrt(k[faceCelli])
                   /(muw[facei]/rhow[facei]);

                if (yPlus > 11.6)
                {
                    mutw[facei] =
                        muw[facei]
                       *(yPlus*kappa_/::log(E_*yPlus) - 1);
                }
                else
                {
                    mutw[facei] = 0.0;
                }
            }
        }
    }
}
