{
    IOobject ioPoints
    (
        "points",
        runTime.timeName(),
        polyMesh::meshSubDir,
        mesh
    );

    if (ioPoints.headerOk())
    {
        // Reading new points
        pointIOField newPoints
        (
            IOobject
            (
                "points",
                mesh.time().timeName(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        mesh.movePoints(newPoints);
    }
}
