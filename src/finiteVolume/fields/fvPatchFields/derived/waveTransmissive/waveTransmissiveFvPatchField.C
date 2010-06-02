/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "waveTransmissiveFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicholsonDdtScheme.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(p, iF),
    psiName_("psi"),
    gamma_(0.0)
{}


template<class Type>
waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const waveTransmissiveFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    advectiveFvPatchField<Type>(ptf, p, iF, mapper),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_)
{}


template<class Type>
waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    advectiveFvPatchField<Type>(p, iF, dict),
    psiName_(dict.lookupOrDefault<word>("psi", "psi")),
    gamma_(readScalar(dict.lookup("gamma")))
{}


template<class Type>
waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const waveTransmissiveFvPatchField& ptpsf
)
:
    advectiveFvPatchField<Type>(ptpsf),
    psiName_(ptpsf.psiName_),
    gamma_(ptpsf.gamma_)
{}


template<class Type>
waveTransmissiveFvPatchField<Type>::waveTransmissiveFvPatchField
(
    const waveTransmissiveFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(ptpsf, iF),
    psiName_(ptpsf.psiName_),
    gamma_(ptpsf.gamma_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<scalarField> waveTransmissiveFvPatchField<Type>::advectionSpeed() const
{
    // Lookup the velocity and compressibility of the patch
    const fvPatchField<scalar>& psip = this->patch().lookupPatchField
    (
        psiName_,
        reinterpret_cast<const volScalarField*>(0),
        reinterpret_cast<const scalar*>(0)
    );

    const surfaceScalarField& phi =
        this->db().objectRegistry::lookupObject<surfaceScalarField>
        (this->phiName_);

    fvsPatchField<scalar> phip = this->patch().lookupPatchField
    (
        this->phiName_,
        reinterpret_cast<const surfaceScalarField*>(0),
        reinterpret_cast<const scalar*>(0)
    );

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchScalarField& rhop = this->patch().lookupPatchField
        (
            this->rhoName_,
            reinterpret_cast<const volScalarField*>(0),
            reinterpret_cast<const scalar*>(0)
        );

        phip /= rhop;
    }

    // Calculate the speed of the field wave w
    // by summing the component of the velocity normal to the boundary
    // and the speed of sound (sqrt(gamma_/psi)).
    return phip/this->patch().magSf() + sqrt(gamma_/psip);
}


template<class Type>
void waveTransmissiveFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    if (this->phiName_ != "phi")
    {
        os.writeKeyword("phi") << this->phiName_ << token::END_STATEMENT << nl;
    }
    if (this->rhoName_ != "rho")
    {
        os.writeKeyword("rho") << this->rhoName_ << token::END_STATEMENT << nl;
    }
    if (psiName_ != "psi")
    {
        os.writeKeyword("psi") << psiName_ << token::END_STATEMENT << nl;
    }
    os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;

    if (this->lInf_ > SMALL)
    {
        os.writeKeyword("fieldInf") << this->fieldInf_
            << token::END_STATEMENT << nl;
        os.writeKeyword("lInf") << this->lInf_
            << token::END_STATEMENT << nl;
    }

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
