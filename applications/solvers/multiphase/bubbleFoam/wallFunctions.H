{
    labelList cellBoundaryFaceCount(epsilon.size(), 0);

    scalar Cmu25 = ::pow(Cmu.value(), 0.25);
    scalar Cmu75 = ::pow(Cmu.value(), 0.75);
    scalar kappa_ = kappa.value();
    scalar nub_ = nub.value();

    const fvPatchList& patches = mesh.boundary();

    //- Initialise the near-wall P field to zero
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                epsilon[faceCelli] = 0.0;
                G[faceCelli] = 0.0;
            }
        }
    }

    //- Accumulate the wall face contributions to epsilon and G
    //  Increment cellBoundaryFaceCount for each face for averaging
    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            const scalarField& nutbw = nutb.boundaryField()[patchi];

            scalarField magFaceGradU(mag(Ub.boundaryField()[patchi].snGrad()));

            forAll(currPatch, facei)
            {
                label faceCelli = currPatch.faceCells()[facei];

                scalar yPlus =
                    Cmu25*y[patchi][facei]
                    *::sqrt(k[faceCelli])
                    /nub_;


                // For corner cells (with two boundary or more faces),
                // epsilon and G in the near-wall cell are calculated
                // as an average

                cellBoundaryFaceCount[faceCelli]++;

                epsilon[faceCelli] +=
                     Cmu75*::pow(k[faceCelli], 1.5)
                    /(kappa_*y[patchi][facei]);

                if (yPlus > 11.6)
                {
                    G[faceCelli] +=
                        (nutbw[facei] + nub_)*magFaceGradU[facei]
                        *Cmu25*::sqrt(k[faceCelli])
                        /(kappa_*y[patchi][facei]);
                }
            }
        }
    }


    // perform the averaging

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label faceCelli = curPatch.faceCells()[facei];

                epsilon[faceCelli] /= cellBoundaryFaceCount[faceCelli];
                G[faceCelli] /= cellBoundaryFaceCount[faceCelli];
            }
        }
    }
}
