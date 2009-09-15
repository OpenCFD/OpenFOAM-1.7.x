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

#include "fanFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(fvPatchScalarField, fanFvPatchScalarField);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Specialisation of the jump-condition for the pressure
template<>
void Foam::fanFvPatchField<Foam::scalar>::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    jump_ = f_[0];

    if (f_.size() > 1)
    {
        const surfaceScalarField& phi =
            db().lookupObject<surfaceScalarField>("phi");

        const fvsPatchField<scalar>& phip =
            patch().patchField<surfaceScalarField, scalar>(phi);

        scalarField Un = max
        (
            scalarField::subField(phip, size()/2)
           /scalarField::subField(patch().magSf(), size()/2),
            scalar(0)
        );

        if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
        {
            Un /=
                scalarField::subField
                (
                    patch().lookupPatchField<volScalarField, scalar>("rho"),
                    size()/2
                );
        }

        for(label i=1; i<f_.size(); i++)
        {
            jump_ += f_[i]*pow(Un, i);
        }

        jump_ = max(jump_, scalar(0));
    }

    jumpCyclicFvPatchField<scalar>::updateCoeffs();
}


// ************************************************************************* //
