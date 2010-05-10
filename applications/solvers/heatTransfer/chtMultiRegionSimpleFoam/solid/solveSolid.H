{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix tEqn
        (
            -fvm::laplacian(K, T)
        );
        tEqn.relax();
        eqnResidual = tEqn.solve().initialResidual();
        maxResidual = max(eqnResidual, maxResidual);

    }

    Info<< "Min/max T:" << min(T).value() << ' '
        << max(T).value() << endl;
}
