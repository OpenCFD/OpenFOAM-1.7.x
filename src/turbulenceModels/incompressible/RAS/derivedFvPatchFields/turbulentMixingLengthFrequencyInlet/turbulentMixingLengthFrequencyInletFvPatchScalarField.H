/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2011 OpenCFD Ltd.
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
    Foam::incompressible::turbulentMixingLengthFrequencyInletFvPatchScalarField

Description
    Calculate omega via the mixing length

    Example of the boundary condition specification:
    @verbatim
        inlet
        {
            type            turbulentMixingLengthFrequencyInlet;
            mixingLength    0.005;         // 5 mm
            k               k;             // turbulent k field
            value           uniform 5;     // initial value
        }
    @endverbatim

SourceFiles
    turbulentMixingLengthFrequencyInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleturbulentMixingLengthFrequencyInletFvPatchScalarField_H
#define incompressibleturbulentMixingLengthFrequencyInletFvPatchScalarField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
    Class turbulentMixingLengthFrequencyInletFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class turbulentMixingLengthFrequencyInletFvPatchScalarField
:
    public inletOutletFvPatchScalarField
{
    // Private data

        //- Turbulent length scale
        scalar mixingLength_;

        //- Name of the flux field
        word phiName_;

        //- Name of the turbulent kinetic energy field
        word kName_;


public:

    //- Runtime type information
    TypeName("turbulentMixingLengthFrequencyInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentMixingLengthFrequencyInletFvPatchScalarField
        //  onto a new patch
        turbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const turbulentMixingLengthFrequencyInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const turbulentMixingLengthFrequencyInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentMixingLengthFrequencyInletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentMixingLengthFrequencyInletFvPatchScalarField
        (
            const turbulentMixingLengthFrequencyInletFvPatchScalarField&,
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
                new turbulentMixingLengthFrequencyInletFvPatchScalarField
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

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
