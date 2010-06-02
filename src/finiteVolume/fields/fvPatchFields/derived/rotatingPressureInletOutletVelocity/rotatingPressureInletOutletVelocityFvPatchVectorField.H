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
    Foam::rotatingPressureInletOutletVelocityFvPatchVectorField

Description
    Velocity inlet/outlet boundary condition in a rotating frame
    for patches where the pressure is specified.  zero-gradient is applied for
    outflow (as defined by the flux) and for inflow the velocity is obtained
    from the flux with a direction normal to the patch faces.

SourceFiles
    rotatingPressureInletOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingPressureInletOutletVelocityFvPatchVectorField_H
#define rotatingPressureInletOutletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "pressureInletOutletVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class rotatingPressureInletOutletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class rotatingPressureInletOutletVelocityFvPatchVectorField
:
    public pressureInletOutletVelocityFvPatchVectorField
{
    // Private data

        //- Angular velocity of the frame
        vector omega_;


    // Private member functions

        //- Calculate the tangentialVelocity from omega
        void calcTangentialVelocity();


public:

    //- Runtime type information
    TypeName("rotatingPressureInletOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        rotatingPressureInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rotatingPressureInletOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  rotatingPressureInletOutletVelocityFvPatchVectorField
        //  onto a new patch
        rotatingPressureInletOutletVelocityFvPatchVectorField
        (
            const rotatingPressureInletOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        rotatingPressureInletOutletVelocityFvPatchVectorField
        (
            const rotatingPressureInletOutletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new rotatingPressureInletOutletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        rotatingPressureInletOutletVelocityFvPatchVectorField
        (
            const rotatingPressureInletOutletVelocityFvPatchVectorField&,
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
                new rotatingPressureInletOutletVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return the angular velocity of rotation
            const vector& omega() const
            {
                return omega_;
            }

            //- Reset the angular velocity of rotation
            //  and update the tangentialVelocity
            void setOmega(const vector& omega)
            {
                omega_ = omega;
                calcTangentialVelocity();
            }


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
