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

Class
    Foam::compressible::RASModels::
    mutSpalartAllmarasWallFunctionFvPatchScalarField

Description
    Wall function boundary condition for walls

SourceFiles
    mutSpalartAllmarasWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleMutSpalartAllmarasWallFunctionFvPatchScalarField_H
#define compressibleMutSpalartAllmarasWallFunctionFvPatchScalarField_H

#include "mutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
      Class mutSpalartAllmarasWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class mutSpalartAllmarasWallFunctionFvPatchScalarField
:
    public mutWallFunctionFvPatchScalarField
{
protected:

    // Protected member functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcMut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField& magGradU) const;


public:

    //- Runtime type information
    TypeName("mutSpalartAllmarasWallFunction");


    // Constructors

        //- Construct from patch and internal field
        mutSpalartAllmarasWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mutSpalartAllmarasWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  mutSpalartAllmarasWallFunctionFvPatchScalarField
        //  onto a new patch
        mutSpalartAllmarasWallFunctionFvPatchScalarField
        (
            const mutSpalartAllmarasWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        mutSpalartAllmarasWallFunctionFvPatchScalarField
        (
            const mutSpalartAllmarasWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new mutSpalartAllmarasWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        mutSpalartAllmarasWallFunctionFvPatchScalarField
        (
            const mutSpalartAllmarasWallFunctionFvPatchScalarField&,
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
                new mutSpalartAllmarasWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
