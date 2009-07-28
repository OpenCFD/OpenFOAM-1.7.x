    volVectorField Ur = Ua - Ub;
    volScalarField magUr = mag(Ur);

    volScalarField CdaMagUr =
        (24.0*nub/da)*(scalar(1) + 0.15*pow(da*magUr/nub, 0.687));

    volScalarField CdbMagUr =
        (24.0*nua/db)*(scalar(1) + 0.15*pow(db*magUr/nua, 0.687));

    volScalarField dragCoef
    (
        "Cd",
        0.75*(beta*rhob*CdaMagUr/da + alpha*rhoa*CdbMagUr/db)
    );

    volVectorField liftCoeff = Cl*(beta*rhob + alpha*rhoa)*(Ur ^ fvc::curl(U));
