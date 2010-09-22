    IOobject phiaHeader
    (
        "phia",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

    autoPtr<surfaceScalarField> phiaPtr(NULL);

    if (phiaHeader.headerOk())
    {
        Info<< "Reading face flux field phia\n" << endl;

        phiaPtr.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phia",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        Info<< "Calculating face flux field phia\n" << endl;

        wordList phiTypes
        (
            Ua.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        for (label i=0; i<Ua.boundaryField().size(); i++)
        {
            if (isA<fixedValueFvPatchVectorField>(Ua.boundaryField()[i]))
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        phiaPtr.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phia",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(Ua) & mesh.Sf(),
                phiTypes
            )
        );
    }

    surfaceScalarField& phia = phiaPtr();
