    Info<< "Create mesh for time = "
        << runTime.timeName() << nl << endl;

    autoPtr<engineMesh> meshPtr
    (
        engineMesh::New
        (
            IOobject
            (
                engineMesh::defaultRegion,
                runTime.timeName(),
                runTime
            )
        )
    );

    engineMesh& mesh = meshPtr();
