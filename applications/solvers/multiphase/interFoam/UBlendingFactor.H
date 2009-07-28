    surfaceScalarField alpha1f = fvc::interpolate(alpha1);
    surfaceScalarField UBlendingFactor
    (
        "UBlendingFactor",
        sqrt(max(min(4*alpha1f*(1.0 - alpha1f), 1.0), 0.0))
    );
