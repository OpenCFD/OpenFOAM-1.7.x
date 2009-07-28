Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

Info<< "Reading inlet velocity  uInfX\n" << endl;

dimensionedScalar uInfX
(
    "uInfx",
    dimensionSet(0, 1, -1, 0, 0),
    U.boundaryField()[3][0].x()
);
Info << "U at inlet = " << uInfX.value() << " m/s" << endl;

dimensionedScalar radius
( 
    "radius",
    dimensionSet(0, 1, 0, 0, 0),
    mag(U.mesh().boundary()[4].Cf()[0])
);
    
Info << "Cylinder radius = " << radius.value() << " m" << endl;

volVectorField UA
(
    IOobject
    (
        "UA",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    U
);
