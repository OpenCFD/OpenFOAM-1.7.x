    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Creating face flux\n" << endl;
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", mesh.Sf().dimensions()*U.dimensions(), 0.0)
    );


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );


    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedVector Ubar
    (
        transportProperties.lookup("Ubar")
    );

    vector flowDirection = (Ubar/mag(Ubar)).value();
    tensor flowMask = sqr(flowDirection);


    // Search for wall patches faces and store normals

    scalar nWallFaces(0);
    vector wallNormal(vector::zero);

    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& currPatch = patches[patchi];

        if (isA<wallFvPatch>(currPatch))
        {
            forAll(currPatch, facei)
            {
                nWallFaces++;

                if (nWallFaces == 1)
                {
                    wallNormal =
                      - mesh.Sf().boundaryField()[patchi][facei]
                        /mesh.magSf().boundaryField()[patchi][facei];
                }
                else if (nWallFaces == 2)
                {
                    vector wallNormal2 =
                        mesh.Sf().boundaryField()[patchi][facei]
                        /mesh.magSf().boundaryField()[patchi][facei];

                    //- Check that wall faces are parallel
                    if
                    (
                        mag(wallNormal & wallNormal2) > 1.01
                      ||mag(wallNormal & wallNormal2) < 0.99
                    )
                    {
                        Info<< "boundaryFoam: wall faces are not parallel"
                            << endl
                            << abort(FatalError);
                    }
                }
                else
                {
                    Info<< "boundaryFoam: number of wall faces > 2"
                        << endl
                        << abort(FatalError);
                }
            }
        }
    }


    //- create position array for graph generation
    scalarField y = wallNormal & mesh.C().internalField();


    dimensionedVector gradP
    (
        "gradP",
        dimensionSet(0, 1, -2, 0, 0),
        vector(0, 0, 0)
    );
