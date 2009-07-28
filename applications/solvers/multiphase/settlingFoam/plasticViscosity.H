volScalarField plasticViscosity
(
    const dimensionedScalar& plasticViscosityCoeff,
    const dimensionedScalar& plasticViscosityExponent,
    const volScalarField& alpha
)
{
    return
    (
         plasticViscosityCoeff*
         (
             pow(10.0, plasticViscosityExponent*alpha + SMALL) - scalar(1)
         )
    );
}
