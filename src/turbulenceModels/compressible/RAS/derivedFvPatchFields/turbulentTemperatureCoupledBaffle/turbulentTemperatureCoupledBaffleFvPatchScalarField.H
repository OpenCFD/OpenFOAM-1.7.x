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

Class
    turbulentTemperatureCoupledBaffleFvPatchScalarField

Description
    Harmonic fixed value boundary condition for temperature, to be used
    for heat-transfer on back-to-back baffles.

    If my temperature is T1, heat conductivity K1 and neighbour is T2,K2

    both sides get fixedValue (K1/dx1*T1 + K2/dx2*T2)/(K1/dx1+K2/dx2)

    Example usage:
        myInterfacePatchName
        {
            type                compressible::turbulentTemperatureCoupledBaffle;
            neighbourFieldName  T;
            K                   K; // or none
            value               uniform 300;
        }

    Needs to be on underlying directMapped(Wall)FvPatch.

    Note: if K is "none" looks up RASModel and basicThermo, otherwise expects
    the solver to calculate a 'K' field.

    Note: runs in parallel with arbitrary decomposition. Uses directMapped
    functionality to calculate exchange.

    Note: lags interface data so both sides use same data.
    - problem: schedule to calculate average would interfere
    with standard processor swaps.
    - so: updateCoeffs sets both to same Twall. Only need to do
    this for last outer iteration but don't have access to this.

SourceFiles
    turbulentTemperatureCoupledBaffleFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentTemperatureCoupledBaffleFvPatchScalarField_H
#define turbulentTemperatureCoupledBaffleFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class turbulentTemperatureCoupledBaffleFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentTemperatureCoupledBaffleFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of field on the neighbour region
        const word neighbourFieldName_;

        //- Name of thermal conductivity field
        const word KName_;


    // Private Member Functions

        //- Am I or neighbour owner of interface
        bool interfaceOwner(const polyMesh&, const polyPatch&) const;


public:

    //- Runtime type information
    TypeName("compressible::turbulentTemperatureCoupledBaffle");


    // Constructors

        //- Construct from patch and internal field
        turbulentTemperatureCoupledBaffleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentTemperatureCoupledBaffleFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleFvPatchScalarField onto a new patch
        turbulentTemperatureCoupledBaffleFvPatchScalarField
        (
            const turbulentTemperatureCoupledBaffleFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentTemperatureCoupledBaffleFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentTemperatureCoupledBaffleFvPatchScalarField
        (
            const turbulentTemperatureCoupledBaffleFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentTemperatureCoupledBaffleFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get corresponding K field
        tmp<scalarField> K() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
