Info<< "magnetic flux divergence error = "
    << runTime.deltaT().value()
      *mag(fvc::div(phiB))().weightedAverage(mesh.V()).value()
    << endl;
