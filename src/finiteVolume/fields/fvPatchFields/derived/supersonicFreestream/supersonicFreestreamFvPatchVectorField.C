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

#include "supersonicFreestreamFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

supersonicFreestreamFvPatchVectorField::supersonicFreestreamFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    UInf_(vector::zero),
    pInf_(0),
    TInf_(0),
    gamma_(0)
{
    refValue() = patchInternalField();
    refGrad() = vector::zero;
    valueFraction() = 1;
}


supersonicFreestreamFvPatchVectorField::supersonicFreestreamFvPatchVectorField
(
    const supersonicFreestreamFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    UInf_(ptf.UInf_),
    pInf_(ptf.pInf_),
    TInf_(ptf.TInf_),
    gamma_(ptf.gamma_)
{}


supersonicFreestreamFvPatchVectorField::supersonicFreestreamFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF),
    UInf_(dict.lookup("UInf")),
    pInf_(readScalar(dict.lookup("pInf"))),
    TInf_(readScalar(dict.lookup("TInf"))),
    gamma_(readScalar(dict.lookup("gamma")))
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<vector>::operator=(patchInternalField());
    }

    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 1;

    if (pInf_ < SMALL)
    {
        FatalIOErrorIn
        (
            "supersonicFreestreamFvPatchVectorField::"
            "supersonicFreestreamFvPatchVectorField"
            "(const fvPatch&, const vectorField&, const dictionary&)",
            dict
        )   << "    unphysical pInf specified (pInf <= 0.0)"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalIOError);
    }
}


supersonicFreestreamFvPatchVectorField::supersonicFreestreamFvPatchVectorField
(
    const supersonicFreestreamFvPatchVectorField& sfspvf
)
:
    mixedFvPatchVectorField(sfspvf),
    UInf_(sfspvf.UInf_),
    pInf_(sfspvf.pInf_),
    TInf_(sfspvf.TInf_),
    gamma_(sfspvf.gamma_)
{}


supersonicFreestreamFvPatchVectorField::supersonicFreestreamFvPatchVectorField
(
    const supersonicFreestreamFvPatchVectorField& sfspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(sfspvf, iF),
    UInf_(sfspvf.UInf_),
    pInf_(sfspvf.pInf_),
    TInf_(sfspvf.TInf_),
    gamma_(sfspvf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void supersonicFreestreamFvPatchVectorField::updateCoeffs()
{
    if (!size() || updated())
    {
        return;
    }

    const fvPatchField<scalar>& pT =
        patch().lookupPatchField<volScalarField, scalar>("T");

    const fvPatchField<scalar>& pp =
        patch().lookupPatchField<volScalarField, scalar>("p");

    const fvPatchField<scalar>& ppsi =
        patch().lookupPatchField<volScalarField, scalar>("psi");

    // Need R of the free-stream flow.  Assume R is independent of location
    // along patch so use face 0
    scalar R = 1.0/(ppsi[0]*pT[0]);

    scalar MachInf = mag(UInf_)/sqrt(gamma_*R*TInf_);

    if (MachInf < 1.0)
    {
        FatalErrorIn
        (
            "supersonicFreestreamFvPatchVectorField::updateCoeffs()"
        )   << "    MachInf < 1.0, free stream must be supersonic"
            << "\n    on patch " << this->patch().name()
            << " of field " << this->dimensionedInternalField().name()
            << " in file " << this->dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    vectorField& Up = refValue();
    valueFraction() = 1;

    // get the near patch internal cell values
    vectorField U = patchInternalField();


    // Find the component of U normal to the free-stream flow and in the
    // plane of the free-stream flow and the patch normal

    // Direction of the free-stream flow
    vector UInfHat = UInf_/mag(UInf_);

    // Normal to the plane defined by the free-stream and the patch normal
    vectorField nnInfHat = UInfHat ^ patch().nf();

    //vectorField UnInf = U - nnInfHat*(nnInfHat & U);
    //vectorField Un = UnInf - UInfHat*(UInfHat & UnInf);
    //vectorField nHatInf =
    //    (Un/(mag(Un) + SMALL)) * sign(patch().nf() & Un);

    // Normal to the free-stream in the plane defined by the free-stream
    // and the patch normal
    vectorField nHatInf = nnInfHat ^ UInfHat;

    // Component of U normal to the free-stream in the plane defined by the
    // free-stream and the patch normal
    vectorField Un = nHatInf*(nHatInf & U);

    // The tangential component is
    vectorField Ut = U - Un;

    // Calculate the Prandtl-Meyer function of the free-stream
    scalar nuMachInf =
        sqrt((gamma_ + 1)/(gamma_ - 1))
       *atan(sqrt((gamma_ - 1)/(gamma_ + 1)*(sqr(MachInf) - 1)))
      - atan(sqr(MachInf) - 1);


    // Set the patch boundary condition based on the Mach number and direction
    // of the flow dictated by the boundary/free-stream pressure difference

    forAll(Up, facei)
    {
        if (pp[facei] >= pInf_) // If outflow
        {
            // Assume supersonic outflow and calculate the boundary velocity
            // according to ???

            scalar fpp =
                sqrt(sqr(MachInf) - 1)
               /(gamma_*sqr(MachInf))*mag(Ut[facei])*log(pp[facei]/pInf_);

            Up[facei] = Ut[facei] + fpp*nHatInf[facei];

            // Calculate the Mach number of the boundary velocity
            scalar Mach = mag(Up[facei])/sqrt(gamma_/ppsi[facei]);

            if (Mach <= 1) // If subsonic
            {
                // Zero-gradient subsonic outflow

                Up[facei] = U[facei];
                valueFraction()[facei] = 0;
            }
        }
        else // if inflow
        {
            // Calculate the Mach number of the boundary velocity
            // from the boundary pressure assuming constant total pressure
            // expansion from the free-stream
            scalar Mach =
                sqrt
                (
                    (2/(gamma_ - 1))*(1 + ((gamma_ - 1)/2)*sqr(MachInf))
                   *pow(pp[facei]/pInf_, (1 - gamma_)/gamma_)
                  - 2/(gamma_ - 1)
                );

            if (Mach > 1) // If supersonic
            {
                // Supersonic inflow is assumed to occur according to the
                // Prandtl-Meyer expansion process

                scalar nuMachf =
                    sqrt((gamma_ + 1)/(gamma_ - 1))
                   *atan(sqrt((gamma_ - 1)/(gamma_ + 1)*(sqr(Mach) - 1)))
                  - atan(sqr(Mach) - 1);

                scalar fpp = (nuMachInf - nuMachf)*mag(Ut[facei]);

                Up[facei] = Ut[facei] + fpp*nHatInf[facei];
            }
            else // If subsonic
            {
                FatalErrorIn
                (
                    "supersonicFreestreamFvPatchVectorField::updateCoeffs()"
                )   << "unphysical subsonic inflow has been generated"
                    << "\n    on patch " << this->patch().name()
                    << " of field " << this->dimensionedInternalField().name()
                    << " in file "
                    << this->dimensionedInternalField().objectPath()
                    << exit(FatalError);
            }
        }
    }

    mixedFvPatchVectorField::updateCoeffs();
}


void supersonicFreestreamFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("UInf") << UInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("pInf") << pInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("TInf") << TInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    supersonicFreestreamFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
