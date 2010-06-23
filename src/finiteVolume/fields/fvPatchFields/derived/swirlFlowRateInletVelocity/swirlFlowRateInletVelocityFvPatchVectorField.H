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
    Foam::swirlFlowRateInletVelocityFvPatchVectorField

Description
    Describes a volumetric/mass flow normal vector boundary condition by its
    magnitude as an integral over its area with a swirl component determined
    by the RPM

    The basis of the patch (volumetric or mass) is determined by the
    dimensions of the flux, phi.
    The current density is used to correct the velocity when applying the
    mass basis.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            swirlFlowRateInletVelocity;
        flowRate        0.2;        // Volumetric/mass flow rate [m3/s or kg/s]
        rpm             100;
    }
    @endverbatim

Note
    - The value is positive inwards

SourceFiles
    swirlFlowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef swirlFlowRateInletVelocityFvPatchVectorField_H
#define swirlFlowRateInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
      Class swirlFlowRateInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class swirlFlowRateInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Inlet integral flow rate
        const scalar flowRate_;

        //- Name of the flux transporting the field
        const word phiName_;

        //- Name of the density field used to normalize the mass flux
        const word rhoName_;

        //- RPM
        const scalar rpm_;


public:

   //- Runtime type information
   TypeName("swirlFlowRateInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  flowRateInletVelocityFvPatchVectorField
        //  onto a new patch
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const swirlFlowRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const swirlFlowRateInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new swirlFlowRateInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        swirlFlowRateInletVelocityFvPatchVectorField
        (
            const swirlFlowRateInletVelocityFvPatchVectorField&,
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
                new swirlFlowRateInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the flux
            scalar flowRate() const
            {
                return flowRate_;
            }


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
