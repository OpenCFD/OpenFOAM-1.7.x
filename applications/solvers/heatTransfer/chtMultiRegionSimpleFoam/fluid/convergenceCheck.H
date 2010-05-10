// check convergence
Info<< "maxResidual: " << maxResidual
    << "  convergence criterion: " << convergenceCriterion
    << endl;

if (maxResidual < convergenceCriterion)
{
    Info<< "reached convergence criterion: " << convergenceCriterion << endl;
    runTime.writeAndEnd();
    Info<< "latestTime = " << runTime.timeName() << endl;
}

