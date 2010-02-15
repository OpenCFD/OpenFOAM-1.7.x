    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<rhoChemistryModel> pChemistry
    (
        rhoChemistryModel::New(mesh)
    );
    rhoChemistryModel& chemistry = pChemistry();

    hsReactionThermo& thermo = chemistry.thermo();

    basicMultiComponentMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    word inertSpecie(thermo.lookup("inertSpecie"));

    if (!composition.contains(inertSpecie))
    {
        FatalErrorIn(args.executable())
            << "Specified inert specie '" << inertSpecie << "' not found in "
            << "species list. Available species:" << composition.species()
            << exit(FatalError);
    }

    volScalarField& p = thermo.p();
    volScalarField& hs = thermo.hs();
    const volScalarField& T = thermo.T();
    const volScalarField& psi = thermo.psi();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    Info<< "\nReading field U\n" << endl;
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

    #include "compressibleCreatePhi.H"

    DimensionedField<scalar, volMesh> kappa
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Creating multi-variate interpolation scheme\n" << endl;
    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    forAll(Y, i)
    {
        fields.add(Y[i]);
    }
    fields.add(hs);

    DimensionedField<scalar, volMesh> chemistrySh
    (
        IOobject
        (
            "chemistry::Sh",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("chemistrySh", dimEnergy/dimTime/dimVolume, 0.0)
    );
