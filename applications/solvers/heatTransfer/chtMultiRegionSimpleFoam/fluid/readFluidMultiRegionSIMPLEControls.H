    dictionary simple = fluidRegions[i].solutionDict().subDict("SIMPLE");

    int nNonOrthCorr = 0;
    if (simple.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(simple.lookup("nNonOrthogonalCorrectors"));
    }

    bool momentumPredictor = true;
    if (simple.found("momentumPredictor"))
    {
        momentumPredictor = Switch(simple.lookup("momentumPredictor"));
    }

    bool fluxGradp = false;
    if (simple.found("fluxGradp"))
    {
        fluxGradp = Switch(simple.lookup("fluxGradp"));
    }

    bool transonic = false;
    if (simple.found("transonic"))
    {
        transonic = Switch(simple.lookup("transonic"));
    }
