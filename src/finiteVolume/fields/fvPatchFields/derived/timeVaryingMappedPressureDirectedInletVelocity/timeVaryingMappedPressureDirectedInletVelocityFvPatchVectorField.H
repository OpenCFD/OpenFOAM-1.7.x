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
    Foam::timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField

Description
    Velocity inlet boundary condition for patches where the pressure is
    specified. The inflow velocity is obtained from the flux with the
    "inletDirection" direction mapped from a time/space-varying table.

SourceFiles
    timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField_H
#define timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField_H

#include "timeVaryingMappedFixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class timeVaryingMappedPressureDirectedInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
:
    public timeVaryingMappedFixedValueFvPatchVectorField
{
    // Private data

        word phiName_;
        word rhoName_;


public:

    //- Runtime type information
    TypeName("timeVaryingMappedPressureDirectedInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
        //  onto a new patch
        timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
        (
            const timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
        (
            const timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new
                timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
        (
            const timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField&,
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
                new
                timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return the name of rho
            const word& rhoName() const
            {
                return rhoName_;
            }

            //- Return reference to the name of rho to allow adjustment
            word& rhoName()
            {
                return rhoName_;
            }

            //- Return the name of phi
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of phi to allow adjustment
            word& phiName()
            {
                return phiName_;
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
                const fvPatchVectorField&,
                const labelList&
            );


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
