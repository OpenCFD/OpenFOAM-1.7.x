const volVectorField::GeometricBoundaryField& Ubf = U.boundaryField();
wordList rhoUboundaryTypes = Ubf.types();

forAll(rhoUboundaryTypes, patchi)
{
    if (Ubf[patchi].fixesValue())
    {
        rhoUboundaryTypes[patchi] = fixedRhoUFvPatchVectorField::typeName;
    }
}
