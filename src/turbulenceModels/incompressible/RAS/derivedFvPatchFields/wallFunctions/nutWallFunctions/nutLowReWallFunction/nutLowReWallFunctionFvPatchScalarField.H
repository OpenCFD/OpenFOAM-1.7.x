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
    Foam::incompressible::RASModels::nutLowReWallFunctionFvPatchScalarField

Description
    Boundary condition for turbulent (kinematic) viscosity when for use with
    low Reynolds number models.

    Sets nut to zero, and provides an access function to calculate y+.

SourceFiles
    nutLowReWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutLowReWallFunctionFvPatchScalarField_H
#define nutLowReWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
           Class nutLowReWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutLowReWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected member functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;


public:

    //- Runtime type information
    TypeName("nutLowReWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutLowReWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutLowReWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutLowReWallFunctionFvPatchScalarField
        //  onto a new patch
        nutLowReWallFunctionFvPatchScalarField
        (
            const nutLowReWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutLowReWallFunctionFvPatchScalarField
        (
            const nutLowReWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutLowReWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutLowReWallFunctionFvPatchScalarField
        (
            const nutLowReWallFunctionFvPatchScalarField&,
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
                new nutLowReWallFunctionFvPatchScalarField(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
