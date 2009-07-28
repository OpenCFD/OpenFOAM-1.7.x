surfaceScalarField::GeometricBoundaryField& phiPatches = phi.boundaryField();
const volVectorField::GeometricBoundaryField& rhoUpatches = rhoU.boundaryField();
const surfaceVectorField::GeometricBoundaryField& SfPatches = mesh.Sf().boundaryField();

forAll(phiPatches, patchI)
{
    if (phi.boundaryField().types()[patchI] == "calculated")
    {
        calculatedFvsPatchScalarField& phiPatch =
            refCast<calculatedFvsPatchScalarField>(phiPatches[patchI]);

        phiPatch == (rhoUpatches[patchI] & SfPatches[patchI]);
    }
}
