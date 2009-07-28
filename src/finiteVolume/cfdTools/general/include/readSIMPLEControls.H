    dictionary simple = mesh.solutionDict().subDict("SIMPLE");

    int nNonOrthCorr =
        simple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    bool momentumPredictor =
        simple.lookupOrDefault<Switch>("momentumPredictor", true);

    bool fluxGradp =
        simple.lookupOrDefault<Switch>("fluxGradp", false);

    bool transonic =
        simple.lookupOrDefault<Switch>("transonic", false);

