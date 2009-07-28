const dictionary& stressControl = mesh.solutionDict().subDict("stressAnalysis");

int nCorr(readInt(stressControl.lookup("nCorrectors")));
scalar convergenceTolerance(readScalar(stressControl.lookup("D")));
Switch compactNormalStress(stressControl.lookup("compactNormalStress"));
