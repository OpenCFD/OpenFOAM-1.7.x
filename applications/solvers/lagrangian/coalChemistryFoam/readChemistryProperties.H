Info<< "Reading chemistry properties\n" << endl;

IOdictionary chemistryProperties
(
    IOobject
    (
        "chemistryProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

Switch turbulentReaction(chemistryProperties.lookup("turbulentReaction"));

dimensionedScalar Cmix("Cmix", dimless, 1.0);

if (turbulentReaction)
{
    chemistryProperties.lookup("Cmix") >> Cmix;
}
