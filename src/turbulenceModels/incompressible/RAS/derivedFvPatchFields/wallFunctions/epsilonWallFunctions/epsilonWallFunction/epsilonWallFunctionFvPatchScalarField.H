/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::incompressible::RASModels::epsilonWallFunctionFvPatchScalarField

Description
    Boundary condition for epsilon when using wall functions
    - calculates epsilon and G
    - epsilon values added directly into the matrix to act as a constraint

SourceFiles
    epsilonWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef epsilonWallFunctionFvPatchScalarField_H
#define epsilonWallFunctionFvPatchScalarField_H

#include "fixedInternalValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
           Class epsilonWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class epsilonWallFunctionFvPatchScalarField
:
    public fixedInternalValueFvPatchField<scalar>
{
    // Private data

        //- Name of velocity field
        word UName_;

        //- Name of turbulence kinetic energy field
        word kName_;

        //- Name of turbulence generation field
        word GName_;

        //- Name of laminar viscosity field
        word nuName_;

        //- Name of turbulent viscosity field
        word nutName_;

        //- Cmu coefficient
        scalar Cmu_;

        //- Von Karman constant
        scalar kappa_;

        //- E coefficient
        scalar E_;


    // Private member functions

        //- Check the type of the patch
        void checkType();


public:

    //- Runtime type information
    TypeName("epsilonWallFunction");


    // Constructors

        //- Construct from patch and internal field
        epsilonWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        epsilonWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  epsilonWallFunctionFvPatchScalarField
        //  onto a new patch
        epsilonWallFunctionFvPatchScalarField
        (
            const epsilonWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        epsilonWallFunctionFvPatchScalarField
        (
            const epsilonWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new epsilonWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        epsilonWallFunctionFvPatchScalarField
        (
            const epsilonWallFunctionFvPatchScalarField&,
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
                new epsilonWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patchField
            virtual void evaluate(const Pstream::commsTypes);


        // I-O

            //- Write
            void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
