Info<< "Creating mass source\n" << endl;
scalarTimeActivatedExplicitSourceList massSource
(
    "mass",
    mesh,
    dimMass/dimTime/dimVolume,
    composition.species()
);


Info<< "Creating momentum source\n" << endl;
vectorTimeActivatedExplicitSourceList momentumSource
(
    "momentum",
    mesh,
    dimMass*dimVelocity/dimTime/dimVolume,
    "U"
);

Info<< "Creating energy source\n" << endl;
scalarTimeActivatedExplicitSourceList energySource
(
    "energy",
    mesh,
    dimEnergy/dimTime/dimVolume,
    "h"
);
