    PtrList<fvMesh> fluidRegions(rp.fluidRegionNames().size());

    forAll(rp.fluidRegionNames(), i)
    {
        Info<< "Create fluid mesh for region " << rp.fluidRegionNames()[i]
            << " for time = " << runTime.timeName() << nl << endl;

        fluidRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    rp.fluidRegionNames()[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
