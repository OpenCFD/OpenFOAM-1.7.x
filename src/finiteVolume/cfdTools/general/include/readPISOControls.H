    dictionary piso = mesh.solutionDict().subDict("PISO");

    int nCorr(readInt(piso.lookup("nCorrectors")));

    int nNonOrthCorr =
        piso.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    bool momentumPredictor =
        piso.lookupOrDefault<Switch>("momentumPredictor", true);

    bool transonic =
        piso.lookupOrDefault<Switch>("transonic", false);

    int nOuterCorr =
        piso.lookupOrDefault<int>("nOuterCorrectors", 1);
