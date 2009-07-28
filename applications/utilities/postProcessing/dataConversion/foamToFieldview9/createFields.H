// Construct List of pointers to all vol fields
int nFields = volScalarNames.size() + 3*volVectorNames.size();
List<volScalarField*> volFieldPtrs
(
    nFields,
    reinterpret_cast<volScalarField*>(0)
);

stringList volFieldNames(nFields);

nFields = 0;
{
    // Load all scalar fields and store ptr to it
    forAll(volScalarNames, fieldI)
    {
        word fieldName = volScalarNames[fieldI];

        // Check if present
        IOobject ioHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (ioHeader.headerOk())
        {
           volFieldPtrs[nFields] = new volScalarField
            (
                ioHeader,
                mesh
            );
        }

        fieldName = getFieldViewName(fieldName);

        volFieldNames[nFields] = fieldName;

        nFields++;
    }


    // Load all  (componenents of) vector fields
    forAll(volVectorNames, fieldI)
    {
        word fieldName = volVectorNames[fieldI];

        // Check if present
        IOobject ioHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (ioHeader.headerOk())
        {
            volVectorField vvf(ioHeader, mesh);

            // X component
            volFieldPtrs[nFields] =
                new volScalarField
                (
                    vvf.component(vector::X)
                );

            // Y component
            volFieldPtrs[nFields+1] =
                new volScalarField
                (
                    vvf.component(vector::Y)
                );

            // Z component
            volFieldPtrs[nFields+2] =
                new volScalarField
                (
                    vvf.component(vector::Z)
                );
        }

        fieldName = getFieldViewName(fieldName);

        volFieldNames[nFields]   = fieldName + ("x;" +  fieldName);
        volFieldNames[nFields+1] = fieldName + "y";
        volFieldNames[nFields+2] = fieldName + "z";

        nFields += 3;
    }
}



//
// Construct List of pointers to all surface fields
//



int nSurfFields = surfScalarNames.size() + 3*surfVectorNames.size();
List<surfaceScalarField*> surfFieldPtrs
(
    nSurfFields,
    reinterpret_cast<surfaceScalarField*>(0)
);

stringList surfFieldNames(nSurfFields);

nSurfFields = 0;
{
    // Load all scalar fields
    forAll(surfScalarNames, fieldI)
    {
        word fieldName = surfScalarNames[fieldI];

        // Check if present
        IOobject ioHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (ioHeader.headerOk())
        {
            surfFieldPtrs[nSurfFields] =
                new surfaceScalarField
                (
                    ioHeader,
                    mesh
                );
        }

        fieldName = getFieldViewName(fieldName);

        surfFieldNames[nSurfFields] = fieldName;

        nSurfFields++;
    }


    // Set (componenents of) vector fields
    forAll(surfVectorNames, fieldI)
    {
        word fieldName = surfVectorNames[fieldI];

        // Check if present
        IOobject ioHeader
        (
            fieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (ioHeader.headerOk())
        {
            surfaceVectorField svf(ioHeader, mesh);

            // X component
            surfFieldPtrs[nSurfFields] =
                new surfaceScalarField
                (
                    svf.component(vector::X)
                );

            // Y component
            surfFieldPtrs[nSurfFields+1] =
                new surfaceScalarField
                (
                    svf.component(vector::Y)
                );

            // Z component
            surfFieldPtrs[nSurfFields+2] =
                new surfaceScalarField
                (
                    svf.component(vector::Z)
                );
        }

        fieldName = getFieldViewName(fieldName);

        surfFieldNames[nSurfFields]   = fieldName + ("x;" + fieldName);
        surfFieldNames[nSurfFields+1] = fieldName + "y";
        surfFieldNames[nSurfFields+2] = fieldName + "z";

        nSurfFields += 3;
    }
}
