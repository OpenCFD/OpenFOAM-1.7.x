    dictionary simple = mesh.solutionDict().subDict("SIMPLE");

    int nNonOrthCorr =
        simple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    bool momentumPredictor =
        simple.lookupOrDefault<Switch>("momentumPredictor", true);

    bool transonic =
        simple.lookupOrDefault<Switch>("transonic", false);
