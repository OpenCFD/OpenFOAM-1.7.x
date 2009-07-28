
    Info<< nl << "Reading field boundaryT" << endl;
    volScalarField boundaryT
    (
        IOobject
        (
            "boundaryT",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field boundaryU" << endl;
    volVectorField boundaryU
    (
        IOobject
        (
            "boundaryU",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field rhoN (number density)" << endl;
    volScalarField rhoN
    (
        IOobject
        (
            "rhoN",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field rhoM (mass density)" << endl;
    volScalarField rhoM
    (
        IOobject
        (
            "rhoM",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field rhoNdsmc (dsmc particle density)" << endl;
    volScalarField dsmcRhoN
    (
        IOobject
        (
            "dsmcRhoN",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field momentum (momentum density)" << endl;
    volVectorField momentum
    (
        IOobject
        (
            "momentum",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field linearKE (linear kinetic energy density)"
        << endl;

    volScalarField linearKE
    (
        IOobject
        (
            "linearKE",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field internalE (internal energy density)" << endl;
    volScalarField internalE
    (
        IOobject
        (
            "internalE",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field iDof (internal degree of freedom density)"
        << endl;

    volScalarField iDof
    (
        IOobject
        (
            "iDof",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field q (surface heat transfer)" << endl;
    volScalarField q
    (
        IOobject
        (
            "q",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Reading field fD (surface force density)" << endl;
    volVectorField fD
    (
        IOobject
        (
            "fD",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< nl << "Constructing dsmcCloud " << endl;

    dsmcCloud dsmc("dsmc", boundaryT, boundaryU);
