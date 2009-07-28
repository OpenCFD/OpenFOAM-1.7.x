Info<< "\nConstructing multi-component mass flow rate point sources" << endl;
timeActivatedExplicitMulticomponentPointSource pointMassSources
(
    "pointMassSources",
    mesh,
    Y,
    dimMass/dimVolume/dimTime
);
