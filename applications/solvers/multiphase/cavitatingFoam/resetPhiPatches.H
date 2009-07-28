fvsPatchScalarFieldField& phiPatches = phi.boundaryField();
const fvPatchScalarFieldField& rhoPatches = rho.boundaryField();
const fvPatchVectorFieldField& Upatches = U.boundaryField();
const fvsPatchVectorFieldField& SfPatches = mesh.Sf().boundaryField();

forAll(phiPatches, patchI)
{
    if (phi.boundaryField().types()[patchI] == "calculated")
    {
        calculatedFvsPatchScalarField& phiPatch =
            refCast<calculatedFvsPatchScalarField>(phiPatches[patchI]);

        phiPatch == ((rhoPatches[patchI]*Upatches[patchI]) & SfPatches[patchI]);
    }
}
