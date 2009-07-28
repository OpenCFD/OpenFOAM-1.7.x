    const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

    int nCorr(readInt(pimple.lookup("nCorrectors")));

    int nNonOrthCorr =
        pimple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    bool momentumPredictor =
        pimple.lookupOrDefault<Switch>("momentumPredictor", true);
