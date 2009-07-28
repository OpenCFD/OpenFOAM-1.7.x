    dimensionedScalar gradP
    (
        "gradP",
        dimensionSet(0, 1, -2, 0, 0),
        0.0
    );


    IFstream gradPFile
    (
        runTime.path()/runTime.timeName()/"uniform"/"gradP.raw"
    );

    if(gradPFile.good())
    {
        gradPFile >> gradP;
        Info<< "Reading average pressure gradient" <<endl
            << endl;
    }
    else
    {
        Info<< "Initializing with 0 pressure gradient" <<endl
            << endl;
    };
