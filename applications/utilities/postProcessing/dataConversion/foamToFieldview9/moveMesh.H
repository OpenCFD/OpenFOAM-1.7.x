
//
// Check if new points (so moving mesh)
//
{
    IOobject pointsHeader
    (
        "points",
        runTime.timeName(),
        polyMesh::defaultRegion,
        runTime
    );
    if (pointsHeader.headerOk())
    {
        // points exists for time step, let's read them
        Info<< "    Points file detected - updating points" << endl;

        // Reading new points
        pointIOField newPoints
        (
            IOobject
            (
                "points",
                runTime.timeName(),
                polyMesh::defaultRegion,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        mesh.polyMesh::movePoints(newPoints);
    }
}
