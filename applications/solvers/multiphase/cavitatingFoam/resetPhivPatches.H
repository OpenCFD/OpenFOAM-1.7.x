surfaceScalarField::GeometricBoundaryField& phivPatches = phiv.boundaryField();
const volVectorField::GeometricBoundaryField& Upatches = U.boundaryField();
const surfaceVectorField::GeometricBoundaryField& SfPatches = mesh.Sf().boundaryField();

forAll(phivPatches, patchI)
{
    if (phiv.boundaryField().types()[patchI] == "calculated")
    {
        calculatedFvsPatchScalarField& phivPatch =
            refCast<calculatedFvsPatchScalarField>(phivPatches[patchI]);

        phivPatch == (Upatches[patchI] & SfPatches[patchI]);
    }
}
