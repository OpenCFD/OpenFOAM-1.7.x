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

Class
    Foam::fixedPressureCompressibleDensityFvPatchScalarField

Description
    Calculate compressible density as a function of pressure and fluid
    properties.

    Example of the boundary condition specification:
    @verbatim
        inlet
        {
            type        fixedPressureCompressibleDensity;
            p           p;               // Name of static pressure field
            value       uniform 1;       // Initial value
        }
    @endverbatim

SourceFiles
    fixedPressureCompressibleDensityFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedPressureCompressibleDensityFvPatchScalarField_H
#define fixedPressureCompressibleDensityFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class fixedPressureCompressibleDensityFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedPressureCompressibleDensityFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of static pressure field
        word pName_;


public:

    //- Runtime type information
    TypeName("fixedPressureCompressibleDensity");


    // Constructors

        //- Construct from patch and internal field
        fixedPressureCompressibleDensityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedPressureCompressibleDensityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  fixedPressureCompressibleDensityFvPatchScalarField
        //  onto a new patch
        fixedPressureCompressibleDensityFvPatchScalarField
        (
            const fixedPressureCompressibleDensityFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedPressureCompressibleDensityFvPatchScalarField
        (
            const fixedPressureCompressibleDensityFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedPressureCompressibleDensityFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedPressureCompressibleDensityFvPatchScalarField
        (
            const fixedPressureCompressibleDensityFvPatchScalarField&,
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
                new fixedPressureCompressibleDensityFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

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
