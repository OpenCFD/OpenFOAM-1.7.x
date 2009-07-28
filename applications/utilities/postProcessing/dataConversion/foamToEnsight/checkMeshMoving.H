// check for "points" in all of the result directories

bool meshMoving = true;
if (Times.size() > 2)
{
    for (label n1=2; n1<Times.size() && meshMoving; ++n1)
    {
        meshMoving = IOobject
        (
            "points",
            Times[n1].name(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ
        ).headerOk();
    }
}
else
{
    meshMoving = false;
}
