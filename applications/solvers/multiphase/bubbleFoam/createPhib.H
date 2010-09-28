    IOobject phibHeader
    (
        "phib",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ
    );

    autoPtr<surfaceScalarField> phibPtr(NULL);

    if (phibHeader.headerOk())
    {
        Info<< "Reading face flux field phib\n" << endl;

        phibPtr.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phib",
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
        Info<< "Calculating face flux field phib\n" << endl;

        wordList phiTypes
        (
            Ub.boundaryField().size(),
            calculatedFvPatchScalarField::typeName
        );

        forAll(Ub.boundaryField(), i)
        {
            if (isA<fixedValueFvPatchVectorField>(Ub.boundaryField()[i]))
            {
                phiTypes[i] = fixedValueFvPatchScalarField::typeName;
            }
        }

        phibPtr.reset
        (
            new surfaceScalarField
            (
                IOobject
                (
                    "phib",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(Ub) & mesh.Sf(),
                phiTypes
            )
        );
    }

    surfaceScalarField& phib = phibPtr();
