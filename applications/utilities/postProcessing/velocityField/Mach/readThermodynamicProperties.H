    Info<< "Reading thermodynamicProperties\n" << endl;

    IOdictionary thermodynamicProperties
    (
        IOobject
        (
            "thermodynamicProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar R
    (
        thermodynamicProperties.lookup("R")
    );

    dimensionedScalar Cv
    (
        thermodynamicProperties.lookup("Cv")
    );
