Info<< nl << "Reading thermophysicalProperties" << endl;

autoPtr<psiChemistryModel> pChemistry
(
    psiChemistryModel::New(mesh)
);
psiChemistryModel& chemistry = pChemistry();

hsCombustionThermo& thermo = chemistry.thermo();

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

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
    ),
    thermo.rho()
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


volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();
volScalarField& hs = thermo.hs();


#include "compressibleCreatePhi.H"

volScalarField kappa
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

Info << "Creating turbulence model.\n" << nl;
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

Info<< "Creating field DpDt\n" << endl;
volScalarField DpDt =
    fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);


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
