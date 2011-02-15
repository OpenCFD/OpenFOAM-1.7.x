/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
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
    Foam::SRFVelocityFvPatchVectorField

Description
    Freestream velocity patch to be used with SRF model
    to apply the appropriate rotation transformation in time and space.

SourceFiles
    SRFVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef SRFFreestreamVelocityFvPatchVectorField_H
#define SRFFreestreamVelocityFvPatchVectorField_H

#include "inletOutletFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class SRFFreestreamVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class SRFFreestreamVelocityFvPatchVectorField
:
    public inletOutletFvPatchVectorField
{
    // Private data

        //- Velocity of the free stream in the absolute frame
        vector UInf_;


public:

    //- Runtime type information
    TypeName("SRFFreestreamVelocity");


    // Constructors

        //- Construct from patch and internal field
        SRFFreestreamVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        SRFFreestreamVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given SRFFreestreamVelocityFvPatchVectorField
        //  onto a new patch
        SRFFreestreamVelocityFvPatchVectorField
        (
            const SRFFreestreamVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        SRFFreestreamVelocityFvPatchVectorField
        (
            const SRFFreestreamVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new SRFFreestreamVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        SRFFreestreamVelocityFvPatchVectorField
        (
            const SRFFreestreamVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new SRFFreestreamVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the velocity at infinity
            const vector& UInf() const
            {
                return UInf_;
            }

            //- Return reference to the velocity at infinity to allow adjustment
            vector& UInf()
            {
                return UInf_;
            }


        // Evaluation functions

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
