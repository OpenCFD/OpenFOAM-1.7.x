    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    volScalarField& p = thermo.p();
    volScalarField& h = thermo.h();
    const volScalarField& psi = thermo.psi();


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

    #include "compressibleCreatePhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

    dimensionedScalar pMin
    (
        mesh.solutionDict().subDict("SIMPLE").lookup("pMin")
    );

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    dimensionedScalar initialMass = fvc::domainIntegrate(rho);

    porousZones pZones(mesh);
    Switch pressureImplicitPorosity(false);

    int nUCorr = 0;
    if (pZones.size())
    {
        // nUCorrectors for pressureImplicitPorosity
        if (mesh.solutionDict().subDict("SIMPLE").found("nUCorrectors"))
        {
            nUCorr = readInt
            (
                mesh.solutionDict().subDict("SIMPLE").lookup("nUCorrectors")
            );
        }

        if (nUCorr > 0)
        {
            pressureImplicitPorosity = true;
        }
    }

