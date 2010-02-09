/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    Foam::compressible::RASModels::mutLowReWallFunctionFvPatchScalarField

Description
    Boundary condition for turbulent (kinematic) viscosity when for use with
    low Reynolds number models.

    Sets mut to zero, and provides an access to calculate y+.

SourceFiles
    mutLowReWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleMutLowReWallFunctionFvPatchScalarField_H
#define compressibleMutLowReWallFunctionFvPatchScalarField_H

#include "mutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
          Class mutLowReWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mutLowReWallFunctionFvPatchScalarField
:
    public mutWallFunctionFvPatchScalarField
{
protected:

    // Private member functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcMut() const;


public:

    //- Runtime type information
    TypeName("mutLowReWallFunction");


    // Constructors

        //- Construct from patch and internal field
        mutLowReWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mutLowReWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mutLowReWallFunctionFvPatchScalarField
        //  onto a new patch
        mutLowReWallFunctionFvPatchScalarField
        (
            const mutLowReWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mutLowReWallFunctionFvPatchScalarField
        (
            const mutLowReWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mutLowReWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mutLowReWallFunctionFvPatchScalarField
        (
            const mutLowReWallFunctionFvPatchScalarField&,
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
                new mutLowReWallFunctionFvPatchScalarField(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
