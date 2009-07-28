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
    Foam::turbulentInletFvPatchField

Description
    Generate a fluctuating inlet condition by adding a random component
    to a reference (mean) field.
    Input:
        referenceField
            Mean field.
        fluctuationScale
            RMS fluctuation, provided as the fraction of the mean field.
        alpha
            temporal correlation factor;
            the fraction of the new random component added to the previous
            time-step (defaults to 0.1).

SourceFiles
    turbulentInletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentInletFvPatchField_H
#define turbulentInletFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class turbulentInletFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class turbulentInletFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private data

        Random ranGen_;
        Type fluctuationScale_;
        Field<Type> referenceField_;
        scalar alpha_;
        label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("turbulentInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentInletFvPatchField
        //  onto a new patch
        turbulentInletFvPatchField
        (
            const turbulentInletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentInletFvPatchField
        (
            const turbulentInletFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new turbulentInletFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentInletFvPatchField
        (
            const turbulentInletFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new turbulentInletFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the fluctuation scale
            const Type& fluctuationScale() const
            {
                return fluctuationScale_;
            }

            //- Return reference to the fluctuation scale to allow adjustment
            Type& fluctuationScale()
            {
                return fluctuationScale_;
            }

            //- Return the reference field
            const Field<Type>& referenceField() const
            {
                return referenceField_;
            }

            //- Return reference to the reference field to allow adjustment
            Field<Type>& referenceField()
            {
                return referenceField_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "turbulentInletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
