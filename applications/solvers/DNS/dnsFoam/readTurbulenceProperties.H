    Info<< "Reading turbulenceProperties\n" << endl;

    IOdictionary turbulenceProperties
    (
        IOobject
        (
            "turbulenceProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    volVectorField force =
        U/dimensionedScalar("dt", dimTime, runTime.deltaT().value());

    Kmesh K(mesh);
    UOprocess forceGen(K, runTime.deltaT().value(), turbulenceProperties);
