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

    dimensionedScalar Cp = Cv + R;

    dimensionedScalar gamma = Cp/Cv;

    dimensionedScalar Pr = dimensionedScalar("Pr", dimless, 1.0);
    if (thermodynamicProperties.found("Pr"))
    {
        Pr = thermodynamicProperties.lookup("Pr");
    }
