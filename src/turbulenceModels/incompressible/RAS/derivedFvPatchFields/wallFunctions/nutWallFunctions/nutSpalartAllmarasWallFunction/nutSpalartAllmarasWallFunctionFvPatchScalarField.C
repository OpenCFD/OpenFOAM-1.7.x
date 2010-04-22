/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "nutSpalartAllmarasWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField>
nutSpalartAllmarasWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const fvPatchVectorField& Uw = rasModel.U().boundaryField()[patchI];
    const scalarField magGradU = mag(Uw.snGrad());
    const scalarField& nuw = rasModel.nu().boundaryField()[patchI];

    return max
    (
        scalar(0),
        sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - nuw
    );
}


tmp<scalarField> nutSpalartAllmarasWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patch().index()];

    const fvPatchVectorField& Uw =
        rasModel.U().boundaryField()[patch().index()];
    const scalarField magUp = mag(Uw.patchInternalField() - Uw);

    const scalarField& nuw = rasModel.nu().boundaryField()[patch().index()];
    const scalarField& nutw = *this;

    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau();

    forAll(uTau, facei)
    {
        scalar magUpara = magUp[facei];

        scalar ut = sqrt((nutw[facei] + nuw[facei])*magGradU[facei]);

        if (ut > VSMALL)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = min(kappa_*magUpara/ut, 50);
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                scalar f =
                    - ut*y[facei]/nuw[facei]
                    + magUpara/ut
                    + 1/E_*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    y[facei]/nuw[facei]
                  + magUpara/sqr(ut)
                  + 1/E_*kUu*fkUu/ut;

                scalar uTauNew = ut + f/df;
                err = mag((ut - uTauNew)/ut);
                ut = uTauNew;

            } while (ut > VSMALL && err > 0.01 && ++iter < 10);
            uTau[facei] = max(0.0, ut);
        }
    }

    return tuTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutSpalartAllmarasWallFunctionFvPatchScalarField::
nutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF)
{}


nutSpalartAllmarasWallFunctionFvPatchScalarField::
nutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const nutSpalartAllmarasWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


nutSpalartAllmarasWallFunctionFvPatchScalarField::
nutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict)
{}


nutSpalartAllmarasWallFunctionFvPatchScalarField::
nutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const nutSpalartAllmarasWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf)
{}


nutSpalartAllmarasWallFunctionFvPatchScalarField::
nutSpalartAllmarasWallFunctionFvPatchScalarField
(
    const nutSpalartAllmarasWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField>
nutSpalartAllmarasWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchI = patch().index();

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patchI];
    const fvPatchVectorField& Uw = rasModel.U().boundaryField()[patchI];
    const scalarField& nuw = rasModel.nu().boundaryField()[patchI];

    return y*calcUTau(mag(Uw.snGrad()))/nuw;
}


void nutSpalartAllmarasWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutSpalartAllmarasWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
