    IOobject phiHeader
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );


    surfaceScalarField* phiPtr;

    if (phiHeader.headerOk())
    {
        Info<< "Reading face flux field phi\n" << endl;

        phiPtr = new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
    }
    else
    {
        Info<< "Calculating face flux field phi\n" << endl;

        wordList phiTypes
        (
            rhoU.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        phiPtr = new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoU) & mesh.Sf(),
            phiTypes
        );
    }

    surfaceScalarField& phi = *phiPtr;
