    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<hhuCombustionThermo> pThermo
    (
        hhuCombustionThermo::New(mesh)
    );
    hhuCombustionThermo& thermo = pThermo();
    basicMultiComponentMixture& composition = thermo.composition();

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

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();
    volScalarField& h = thermo.h();
    volScalarField& hu = thermo.hu();

    volScalarField& b = composition.Y("b");
    Info<< "min(b) = " << min(b).value() << endl;

    //const volScalarField& T = thermo->T();


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

#   include "compressibleCreatePhi.H"

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

    Info<< "Creating field DpDt\n" << endl;
    volScalarField DpDt
    (
        "DpDt",
        fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p)
    );


    Info<< "Creating the unstrained laminar flame speed\n" << endl;
    autoPtr<laminarFlameSpeed> unstrainedLaminarFlameSpeed
    (
        laminarFlameSpeed::New(thermo)
    );


    Info<< "Reading strained laminar flame speed field Su\n" << endl;
    volScalarField Su
    (
        IOobject
        (
            "Su",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field betav\n" << endl;
    volScalarField betav
    (
        IOobject
        (
            "betav",
            runTime.findInstance(polyMesh::meshSubDir, "betav"),
            polyMesh::meshSubDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    IOdictionary PDRProperties
    (
        IOobject
        (
            "PDRProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    //- Create the drag model
    autoPtr<PDRDragModel> drag = PDRDragModel::New
    (
        PDRProperties,
        turbulence,
        rho,
        U,
        phi
    );

    //- Create the flame-wrinkling model
    autoPtr<XiModel> flameWrinkling = XiModel::New
    (
        PDRProperties,
        thermo,
        turbulence,
        Su,
        rho,
        b,
        phi
    );

    Info<< "Calculating turbulent flame speed field St\n" << endl;
    volScalarField St
    (
        IOobject
        (
            "St",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        flameWrinkling->Xi()*Su
    );


    multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

    if (composition.contains("ft"))
    {
        fields.add(composition.Y("ft"));
    }

    fields.add(b);
    fields.add(h);
    fields.add(hu);
    flameWrinkling->addXi(fields);
