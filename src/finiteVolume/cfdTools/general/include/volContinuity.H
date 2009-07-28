{
    volScalarField conserve = -fvc::div(mesh.phi());
    // The ddt term constructed by hand because it would be wrong for
    // Backward Differencing in time.

    conserve.internalField() +=
        (1.0 - mesh.V0()/mesh.V())/runTime.deltaT().value();

    scalar sumLocalContErr = runTime.deltaT().value()*
        mag(conserve)().weightedAverage(mesh.V()).value();

    scalar globalContErr = runTime.deltaT().value()*
        conserve.weightedAverage(mesh.V()).value();

    Info<< "volume continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr << endl;
}
