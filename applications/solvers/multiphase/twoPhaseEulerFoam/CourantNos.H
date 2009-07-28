#   include "CourantNo.H"

{
    scalar UrCoNum = max
    (
        mesh.surfaceInterpolation::deltaCoeffs()*mag(phia - phib)
       /mesh.magSf()
    ).value()*runTime.deltaT().value();

    Info<< "Max Ur Courant Number = " << UrCoNum << endl;

    CoNum = max(CoNum, UrCoNum);
}
