/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "maxwellSlipUFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFixedValueSlipFvPatchVectorField(p, iF),
    accommodationCoeff_(1.0),
    Uwall_(p.size(), vector(0.0, 0.0, 0.0)),
    thermalCreep_(true),
    curvature_(true)
{}


maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const maxwellSlipUFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFixedValueSlipFvPatchVectorField(tdpvf, p, iF, mapper),
    accommodationCoeff_(tdpvf.accommodationCoeff_),
    Uwall_(tdpvf.Uwall_),
    thermalCreep_(tdpvf.thermalCreep_),
    curvature_(tdpvf.curvature_)
{}


maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFixedValueSlipFvPatchVectorField(p, iF),
    accommodationCoeff_(readScalar(dict.lookup("accommodationCoeff"))),
    Uwall_("Uwall", dict, p.size()),
    thermalCreep_(dict.lookupOrDefault("thermalCreep", true)),
    curvature_(dict.lookupOrDefault("curvature", true))
{
    if
    (
        mag(accommodationCoeff_) < SMALL
        ||
        mag(accommodationCoeff_) > 2.0
    )
    {
        FatalIOErrorIn
        (
            "maxwellSlipUFvPatchScalarField::"
            "maxwellSlipUFvPatchScalarField"
            "(const fvPatch&, const scalarField&, const dictionary&)",
            dict
        )   << "unphysical accommodationCoeff_ specified"
            << "(0 < accommodationCoeff_ <= 1)" << endl
            << exit(FatalError);
    }

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        mixedFixedValueSlipFvPatchVectorField::evaluate();
    }
}


maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const maxwellSlipUFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFixedValueSlipFvPatchVectorField(tdpvf, iF),
    accommodationCoeff_(tdpvf.accommodationCoeff_),
    Uwall_(tdpvf.Uwall_),
    thermalCreep_(tdpvf.thermalCreep_),
    curvature_(tdpvf.curvature_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update the coefficients associated with the patch field
void maxwellSlipUFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchScalarField& pmu =
        patch().lookupPatchField<volScalarField, scalar>("mu");
    const fvPatchScalarField& prho =
        patch().lookupPatchField<volScalarField, scalar>("rho");
    const fvPatchField<scalar>& ppsi =
        patch().lookupPatchField<volScalarField, scalar>("psi");

    Field<scalar> C1 = sqrt(ppsi*mathematicalConstant::pi/2.0)
        *(2.0 - accommodationCoeff_)/accommodationCoeff_;

    Field<scalar> pnu = pmu/prho;
    valueFraction() = (1.0/(1.0 + patch().deltaCoeffs()*C1*pnu));

    refValue() = Uwall_;

    if(thermalCreep_)
    {
        const volScalarField& vsfT =
            this->db().objectRegistry::lookupObject<volScalarField>("T");
        label patchi = this->patch().index();
        const fvPatchScalarField& pT = vsfT.boundaryField()[patchi];
        Field<vector> gradpT = fvc::grad(vsfT)().boundaryField()[patchi];
        vectorField n = patch().nf();

        refValue() -= 3.0*pnu/(4.0*pT)*transform(I - n*n, gradpT);
    }

    if(curvature_)
    {
        const fvPatchTensorField& ptauMC =
            patch().lookupPatchField<volTensorField, tensor>("tauMC");
        vectorField n = patch().nf();

        refValue() -= C1/prho*transform(I - n*n, (n & ptauMC));
    }

    mixedFixedValueSlipFvPatchVectorField::updateCoeffs();
}


// Write
void maxwellSlipUFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("accommodationCoeff")
        << accommodationCoeff_ << token::END_STATEMENT << nl;
    Uwall_.writeEntry("Uwall", os);
    os.writeKeyword("thermalCreep")
        << thermalCreep_ << token::END_STATEMENT << nl;
    os.writeKeyword("curvature") << curvature_ << token::END_STATEMENT << nl;

    os.writeKeyword("refValue")
        << refValue() << token::END_STATEMENT << nl;
    os.writeKeyword("valueFraction")
        << valueFraction() << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, maxwellSlipUFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
