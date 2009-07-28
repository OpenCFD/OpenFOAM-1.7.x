Info<< "\nConstructing coal cloud" << endl;
CoalCloud<gasThermoPhysics> coalParcels
(
    "coalCloud1",
    rho,
    U,
    g,
    thermo
);

Info<< "\nConstructing limestone cloud" << endl;
basicThermoCloud limestoneParcels
(
    "limestoneCloud1",
    rho,
    U,
    g,
    thermo
);
