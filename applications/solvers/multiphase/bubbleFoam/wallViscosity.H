{
    scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    scalar kappa_ = kappa.value();
    scalar E_ = E.value();
    scalar nub_ = nub.value();

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            scalarField& nutw = nutb.boundaryField()[patchi];

            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                // calculate yPlus
                scalar yPlus =
                    Cmu25*y[patchi][facei]
                   *::sqrt(k[faceCelli])
                   /nub_;

                if (yPlus > 11.6)
                {
                    nutw[facei] =
                        yPlus*nub_*kappa_
                       /::log(E_*yPlus)
                      - nub_;
                }
                else
                {
                    nutw[facei] = 0.0;
                }
            }
        }
    }
}
