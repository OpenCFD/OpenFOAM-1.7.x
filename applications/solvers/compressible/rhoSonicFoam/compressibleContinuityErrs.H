{
#   include "rhoEqn.H"
}
{
    scalar sumLocalContErr = (sum(mag(rho - psi*p))/sum(rho)).value();
    scalar globalContErr = (sum(rho - psi*p)/sum(rho)).value();
    cumulativeContErr += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
         << ", global = " << globalContErr
         << ", cumulative = " << cumulativeContErr << endl;
}
