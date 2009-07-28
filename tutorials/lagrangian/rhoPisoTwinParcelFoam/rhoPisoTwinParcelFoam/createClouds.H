    Info<< "Constructing thermoCloud1" << endl;
    basicThermoCloud thermoCloud1
    (
        "thermoCloud1",
        rho,
        U,
        g,
        thermo
    );

    Info<< "Constructing kinematicCloud1" << endl;
    basicKinematicCloud kinematicCloud1
    (
        "kinematicCloud1",
        rho,
        U,
        thermo.mu(),
        g
    );

