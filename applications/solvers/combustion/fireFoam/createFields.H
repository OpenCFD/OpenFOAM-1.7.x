Info<< "Reading thermophysical properties\n" << endl;

autoPtr<hsCombustionThermo> pThermo
(
    hsCombustionThermo::New(mesh)
);

hsCombustionThermo& thermo = pThermo();

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

dimensionedScalar stoicRatio
(
    thermo.lookup("stoichiometricAirFuelMassRatio")
);

volScalarField& p = thermo.p();
volScalarField& hs = thermo.hs();

const volScalarField& psi = thermo.psi();

volScalarField& ft = composition.Y("ft");
volScalarField& fu = composition.Y("fu");

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

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New(rho, U, phi, thermo)
);

IOdictionary combustionProperties
(
    IOobject
    (
        "combustionProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

Info<< "Creating combustion model\n" << endl;
autoPtr<combustionModel> combustion
(
    combustionModel::combustionModel::New
    (
        combustionProperties,
        thermo,
        turbulence(),
        phi,
        rho
    )
);

volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimMass/pow3(dimTime)/dimLength, 0.0)
);


Info<< "Creating field DpDt\n" << endl;
volScalarField DpDt =
    fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);


Info<< "Calculating field g.h\n" << endl;
volScalarField gh("gh", g & mesh.C());

surfaceScalarField ghf("gh", g & mesh.Cf());

p += rho*gh;

thermo.correct();

dimensionedScalar initialMass = fvc::domainIntegrate(rho);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

if (composition.contains("ft"))
{
    fields.add(composition.Y("ft"));
}

if (composition.contains("fu"))
{
    fields.add(composition.Y("fu"));
}

fields.add(hs);
