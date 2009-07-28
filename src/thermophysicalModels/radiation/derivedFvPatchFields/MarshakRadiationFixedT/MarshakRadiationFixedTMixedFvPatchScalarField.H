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
    Foam::MarshakRadiationFixedTMixedFvPatchScalarField

Description
    Marshak boundary for radiation G field
    - radiation temperature specified

SourceFiles
    MarshakRadiationFixedTMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef MarshakRadiationFixedTMixedFvPatchField_H
#define MarshakRadiationFixedTMixedFvPatchField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class MarshakRadiationFixedTMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class MarshakRadiationFixedTMixedFvPatchScalarField
:
    public mixedFvPatchScalarField
{

    // Private data

        //- Radiation temperature field
        scalarField Trad_;

        //- Emissivity
        scalar emissivity_;


public:

    //- Runtime type information
    TypeName("MarshakRadiationFixedT");


    // Constructors

        //- Construct from patch and internal field
        MarshakRadiationFixedTMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        MarshakRadiationFixedTMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given MarshakRadiationFvPatchField onto a new
        //  patch
        MarshakRadiationFixedTMixedFvPatchScalarField
        (
            const MarshakRadiationFixedTMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        MarshakRadiationFixedTMixedFvPatchScalarField
        (
            const MarshakRadiationFixedTMixedFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new MarshakRadiationFixedTMixedFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        MarshakRadiationFixedTMixedFvPatchScalarField
        (
            const MarshakRadiationFixedTMixedFvPatchScalarField&,
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
                new MarshakRadiationFixedTMixedFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the radiation temperature
            const scalarField& Trad() const
            {
                return Trad_;
            }

            //- Return reference to the radiation temperature to allow
            //  adjustment
            scalarField& Trad()
            {
                return Trad_;
            }

            //- Return the emissivity
            const scalar& emissivity() const
            {
                return emissivity_;
            }

            //- Return reference to the emissivity to allow adjustment
            scalar& emissivity()
            {
                return emissivity_;
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
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
