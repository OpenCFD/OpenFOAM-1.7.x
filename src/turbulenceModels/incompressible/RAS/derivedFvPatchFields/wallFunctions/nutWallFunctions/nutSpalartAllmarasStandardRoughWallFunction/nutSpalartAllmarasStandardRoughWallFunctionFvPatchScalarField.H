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
    Foam::incompressible::RASModels::
    nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField

Description
    Wall function boundary condition for rough walls

SourceFiles
    nutSpalartAllmarasStandardWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField_H
#define nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
Class nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
    // Private data

        // Roughness model parameters

            //- Height
            scalar roughnessHeight_;

            //- Constant
            scalar roughnessConstant_;

            //- Scale factor
            scalar roughnessFudgeFactor_;


    // Protected member functions

        //- Calculate yPLus
        virtual tmp<scalarField> calcYPlus(const scalarField& magUp) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;


public:

    //- Runtime type information
    TypeName("nutSpalartAllmarasStandardRoughWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
        //  onto a new patch
        nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
        (
            const nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
        (
            const nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
        (
            const nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField&,
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
                new nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return the roughness height
            const scalar& roughnessHeight() const
            {
                return roughnessHeight_;
            }

            //- Return reference to the roughness height to allow adjustment
            scalar& roughnessHeight()
            {
                return roughnessHeight_;
            }


            //- Return the roughness constant scale
            const scalar& roughnessConstant() const
            {
                return roughnessConstant_;
            }

            //- Return reference to the roughness constant to allow adjustment
            scalar& roughnessConstant()
            {
                return roughnessConstant_;
            }

            //- Return the roughness scale factor
            const scalar& roughnessFudgeFactor() const
            {
                return roughnessFudgeFactor_;
            }

            //- Return reference to the roughness scale factor to allow
            //  adjustment
            scalar& roughnessFudgeFactor()
            {
                return roughnessFudgeFactor_;
            }


    // I-O

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
