    const fvMesh& mesh = fluidRegions[i];

    basicPsiThermo& thermo = thermoFluid[i];
    volScalarField& rho = rhoFluid[i];
    volScalarField& K = KFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    const dimensionedVector& g = gFluid[i];

    compressible::turbulenceModel& turb = turbulence[i];
    volScalarField& DpDt = DpDtFluid[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();
    volScalarField& h = thermo.h();

    const dimensionedScalar massIni("massIni", dimMass, initialMassFluid[i]);
