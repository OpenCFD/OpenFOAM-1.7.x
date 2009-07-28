Info << "Constructing Spray" << endl;

PtrList<gasThermoPhysics> gasProperties(Y.size());
forAll(gasProperties, i)
{
    gasProperties.set
    (
        i,
        new gasThermoPhysics
        (
            dynamic_cast<const reactingMixture<gasThermoPhysics>&>
                (thermo).speciesData()[i]
        )
    );
}

spray dieselSpray
(
    U,
    rho,
    p,
    T,
    composition,
    gasProperties,
    thermo,
    g
);

scalar gasMass0 = fvc::domainIntegrate(rho).value();

if (dieselSpray.twoD())
{
    gasMass0 *= 2.0*mathematicalConstant::pi/dieselSpray.angleOfWedge();
}

gasMass0 -=
    dieselSpray.injectedMass(runTime.value()) - dieselSpray.liquidMass();
