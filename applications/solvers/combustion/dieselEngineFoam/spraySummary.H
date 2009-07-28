    label Nparcels = dieselSpray.size();
    reduce(Nparcels, sumOp<label>());

    Info<< "\nNumber of parcels in system.... | "
        << Nparcels << endl
        << "Injected liquid mass........... | "
        << 1e6*dieselSpray.injectedMass(runTime.value()) << " mg" << endl
        << "Liquid Mass in system.......... | "
        << 1e6*dieselSpray.liquidMass() << " mg" << endl
        << "SMD, Dmax...................... | "
        << dieselSpray.smd()*1e6 << " mu, "
        << dieselSpray.maxD()*1e6 << " mu"
        << endl;

    scalar evapMass =
        dieselSpray.injectedMass(runTime.value())
    - dieselSpray.liquidMass();

    scalar gasMass = fvc::domainIntegrate(rho).value();

    if (dieselSpray.twoD())
    {
        gasMass *= 2.0*mathematicalConstant::pi/dieselSpray.angleOfWedge();
    }

    scalar addedMass = gasMass - gasMass0;

    Info<< "Added gas mass................. | " << 1e6*addedMass << " mg"
        << nl << "Evaporation Continuity Error... | "
        << 1e6*(addedMass - evapMass) << " mg" << endl;
