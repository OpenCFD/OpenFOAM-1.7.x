Info<< "Reading thermophysicalProperties\n" << endl;

// Pr defined as a separate constant to enable calculation of k, currently
// inaccessible through thermo
IOdictionary thermophysicalProperties
(
    IOobject
    (
        "thermophysicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar Pr = dimensionedScalar("Pr", dimless, 1.0);

if (thermophysicalProperties.found("Pr"))
{
    Pr = thermophysicalProperties.lookup("Pr");
}
