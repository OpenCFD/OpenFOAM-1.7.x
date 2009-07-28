fvVectorMatrix UaEqn(Ua, Ua.dimensions()*dimVol/dimTime);
fvVectorMatrix UbEqn(Ub, Ub.dimensions()*dimVol/dimTime);

{
    {
        volTensorField gradUaT = fvc::grad(Ua)().T();
        volTensorField Rca
        (
            "Rca",
            ((2.0/3.0)*I)*(sqr(Ct)*k + nuEffa*tr(gradUaT)) - nuEffa*gradUaT
        );

        if (kineticTheory.on())
        {
            Rca -= ((kineticTheory.lambda()/rhoa)*tr(gradUaT))*tensor(I);
        }

        surfaceScalarField phiRa =
            -fvc::interpolate(nuEffa)*mesh.magSf()*fvc::snGrad(alpha)
            /fvc::interpolate(alpha + scalar(0.001));

        UaEqn =
        (
            (scalar(1) + Cvm*rhob*beta/rhoa)*
            (
                fvm::ddt(Ua)
              + fvm::div(phia, Ua, "div(phia,Ua)")
              - fvm::Sp(fvc::div(phia), Ua)
            )

          - fvm::laplacian(nuEffa, Ua)
          + fvc::div(Rca)

          + fvm::div(phiRa, Ua, "div(phia,Ua)")
          - fvm::Sp(fvc::div(phiRa), Ua)
          + (fvc::grad(alpha)/(fvc::average(alpha) + scalar(0.001)) & Rca)
         ==
        //  g                          // Buoyancy term transfered to p-equation
          - fvm::Sp(beta/rhoa*K, Ua)
        //+ beta/rhoa*K*Ub             // Explicit drag transfered to p-equation
          - beta/rhoa*(liftCoeff - Cvm*rhob*DDtUb)
        );

        UaEqn.relax();
    }

    {
        volTensorField gradUbT = fvc::grad(Ub)().T();
        volTensorField Rcb
        (
            "Rcb",
            ((2.0/3.0)*I)*(k + nuEffb*tr(gradUbT)) - nuEffb*gradUbT
        );

        surfaceScalarField phiRb =
            -fvc::interpolate(nuEffb)*mesh.magSf()*fvc::snGrad(beta)
            /fvc::interpolate(beta + scalar(0.001));

        UbEqn = 
        (
            (scalar(1) + Cvm*rhob*alpha/rhob)*
            (
                fvm::ddt(Ub)
              + fvm::div(phib, Ub, "div(phib,Ub)")
              - fvm::Sp(fvc::div(phib), Ub)
            )

          - fvm::laplacian(nuEffb, Ub)
          + fvc::div(Rcb)

          + fvm::div(phiRb, Ub, "div(phib,Ub)")
          - fvm::Sp(fvc::div(phiRb), Ub)

          + (fvc::grad(beta)/(fvc::average(beta) + scalar(0.001)) & Rcb)
         ==
        //  g                          // Buoyancy term transfered to p-equation
          - fvm::Sp(alpha/rhob*K, Ub)
        //+ alpha/rhob*K*Ua            // Explicit drag transfered to p-equation
          + alpha/rhob*(liftCoeff + Cvm*rhob*DDtUa)
        );

        UbEqn.relax();
    }
}
