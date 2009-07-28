const volScalarField::GeometricBoundaryField& Tbf = T.boundaryField();
wordList rhoEboundaryTypes = Tbf.types();

forAll(rhoEboundaryTypes, patchi)
{
    if (Tbf[patchi].fixesValue())
    {
        rhoEboundaryTypes[patchi] = fixedRhoEFvPatchScalarField::typeName;
    }
    else
    {
        rhoEboundaryTypes[patchi] = mixedRhoEFvPatchScalarField::typeName;
    }
}
