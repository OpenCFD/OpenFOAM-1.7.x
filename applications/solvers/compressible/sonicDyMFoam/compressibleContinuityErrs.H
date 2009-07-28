{
    #include "rhoEqn.H"
}
{
    scalar sumLocalContErr =
        sum
        (
            mag(rho.internalField() - (psi*p)().internalField())
        )/sum(rho.internalField());

    scalar globalContErr =
        sum(rho.internalField() - (psi*p)().internalField())
       /sum(rho.internalField());

    cumulativeContErr += globalContErr;

    Info<< "time step continuity errors : sum local = " << sumLocalContErr
         << ", global = " << globalContErr
         << ", cumulative = " << cumulativeContErr << endl;
}
