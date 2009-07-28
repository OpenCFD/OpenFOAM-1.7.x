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

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Calculating field alpha1\n" << endl;
    volScalarField alpha2("alpha2", scalar(1) - alpha1);

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

    #include "createPhi.H"


    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi);

    dimensionedScalar rho10
    (
        twoPhaseProperties.subDict
        (
            twoPhaseProperties.phase1Name()
        ).lookup("rho0")
    );

    dimensionedScalar rho20
    (
        twoPhaseProperties.subDict
        (
            twoPhaseProperties.phase2Name()
        ).lookup("rho0")
    );

    dimensionedScalar psi1
    (
        twoPhaseProperties.subDict
        (
            twoPhaseProperties.phase1Name()
        ).lookup("psi")
    );

    dimensionedScalar psi2
    (
        twoPhaseProperties.subDict
        (
            twoPhaseProperties.phase2Name()
        ).lookup("psi")
    );

    dimensionedScalar pMin(twoPhaseProperties.lookup("pMin"));

    volScalarField rho1 = rho10 + psi1*p;
    volScalarField rho2 = rho20 + psi2*p;

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
        alpha1*rho1 + alpha2*rho2
    );


    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phi
    );

    volScalarField dgdt =
        pos(alpha2)*fvc::div(phi)/max(alpha2, scalar(0.0001));

    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);

    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );


    wordList pcorrTypes
    (
        p.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );

    for (label i=0; i<p.boundaryField().size(); i++)
    {
        if (p.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }
