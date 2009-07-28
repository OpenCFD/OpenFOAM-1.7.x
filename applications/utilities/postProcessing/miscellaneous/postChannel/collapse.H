    scalarField UMeanXvalues = channelIndexing.collapse
    (
        UMean.component(vector::X)()
    );

    scalarField UMeanYvalues = channelIndexing.collapse
    (
        UMean.component(vector::Y)()
    );

    scalarField UMeanZvalues = channelIndexing.collapse
    (
        UMean.component(vector::Z)()
    );

    scalarField RxxValues = channelIndexing.collapse(Rxx);
    scalarField RyyValues = channelIndexing.collapse(Ryy);
    scalarField RzzValues = channelIndexing.collapse(Rzz);
    scalarField RxyValues = channelIndexing.collapse(Rxy, true);

    scalarField pPrime2MeanValues = channelIndexing.collapse(pPrime2Mean);

    /*
    scalarField epsilonValues = channelIndexing.collapse(epsilonMean);

    scalarField nuMeanValues = channelIndexing.collapse(nuMean);
    scalarField nuPrimeValues = channelIndexing.collapse(nuPrime);

    scalarField gammaDotMeanValues = channelIndexing.collapse(gammaDotMean);
    scalarField gammaDotPrimeValues = channelIndexing.collapse(gammaDotPrime);
    */

    scalarField urmsValues = sqrt(mag(RxxValues));
    scalarField vrmsValues = sqrt(mag(RyyValues));
    scalarField wrmsValues = sqrt(mag(RzzValues));

    scalarField kValues =
        0.5*(sqr(urmsValues) + sqr(vrmsValues) + sqr(wrmsValues));


    const scalarField& y = channelIndexing.y();

    makeGraph(y, UMeanXvalues, "Uf", UMean.path(), gFormat);
    makeGraph(y, urmsValues, "u", UMean.path(), gFormat);
    makeGraph(y, vrmsValues, "v", UMean.path(), gFormat);
    makeGraph(y, wrmsValues, "w", UMean.path(), gFormat);
    makeGraph(y, RxyValues, "uv", UMean.path(), gFormat);
    makeGraph(y, kValues, "k", UMean.path(), gFormat);

    makeGraph(y, pPrime2MeanValues, "pPrime2Mean", UMean.path(), gFormat);

    /*
    makeGraph(y, epsilonValues, "epsilon", UMean.path(), gFormat);
    makeGraph(y, nuMeanValues, "nu", UMean.path(), gFormat);
    makeGraph(y, nuPrimeValues, "nuPrime", UMean.path(), gFormat);
    makeGraph(y, gammaDotMeanValues, "gammaDot", UMean.path(), gFormat);
    makeGraph(y, gammaDotPrimeValues, "gammaDotPrime", UMean.path(), gFormat);
    */
