    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField gamma
    (
        IOobject
        (
            "gamma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        max(min((rho - rholSat)/(rhovSat - rholSat), scalar(1)), scalar(0))
    );
    gamma.oldTime();

    Info<< "Creating compressibilityModel\n" << endl;
    autoPtr<barotropicCompressibilityModel> psiModel =
        barotropicCompressibilityModel::New
        (
            thermodynamicProperties,
            gamma
        );

    const volScalarField& psi = psiModel->psi();

    rho == max
    (
        psi*p
      + (1.0 - gamma)*rhol0
      + ((gamma*psiv + (1.0 - gamma)*psil) - psi)*pSat,
        rhoMin
    );

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

    #include "createPhiv.H"
    #include "compressibleCreatePhi.H"

    Info<< "Reading transportProperties\n" << endl;

    twoPhaseMixture twoPhaseProperties(U, phiv, "gamma");

    // Create incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phiv, twoPhaseProperties)
    );
