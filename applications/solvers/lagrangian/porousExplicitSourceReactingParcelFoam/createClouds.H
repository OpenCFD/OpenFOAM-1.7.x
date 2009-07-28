Info<< "\nConstructing reacting cloud" << endl;
BasicReactingCloud<icoPoly8ThermoPhysics> parcels
(
    "reactingCloud1",
    rho,
    U,
    g,
    thermo
);
