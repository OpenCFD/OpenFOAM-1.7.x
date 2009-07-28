    dictionary piso = fluidRegions[i].solutionDict().subDict("PISO");

    int nCorr(readInt(piso.lookup("nCorrectors")));

    int nNonOrthCorr = 0;
    if (piso.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(piso.lookup("nNonOrthogonalCorrectors"));
    }

    bool momentumPredictor = true;
    if (piso.found("momentumPredictor"))
    {
        momentumPredictor = Switch(piso.lookup("momentumPredictor"));
    }

    bool transonic = false;
    if (piso.found("transonic"))
    {
        transonic = Switch(piso.lookup("transonic"));
    }

    int nOuterCorr = 1;
    if (piso.found("nOuterCorrectors"))
    {
        nOuterCorr = readInt(piso.lookup("nOuterCorrectors"));
    }
