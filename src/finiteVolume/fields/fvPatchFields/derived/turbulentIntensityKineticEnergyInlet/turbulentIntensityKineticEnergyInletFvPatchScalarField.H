/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2009 OpenCFD Ltd.
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
    Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField

Description
    Calculate turbulent kinetic energy from the intensity provided as a
    fraction of the mean velocity

    Example of the boundary condition specification:
    @verbatim
        inlet
        {
            type        turbulentIntensityKineticEnergyInlet;
            intensity   0.05;           // 5% turbulence
            value       uniform 1;      // placeholder
        }
    @endverbatim

SourceFiles
    turbulentIntensityKineticEnergyInletFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentIntensityKineticEnergyInletFvPatchScalarField_H
#define turbulentIntensityKineticEnergyInletFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class turbulentIntensityKineticEnergyInletFvPatch Declaration
\*---------------------------------------------------------------------------*/

class turbulentIntensityKineticEnergyInletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Turbulent intensity as fraction of mean velocity
        scalar intensity_;

public:

    //- Runtime type information
    TypeName("turbulentIntensityKineticEnergyInlet");


    // Constructors

        //- Construct from patch and internal field
        turbulentIntensityKineticEnergyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentIntensityKineticEnergyInletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentIntensityKineticEnergyInletFvPatchScalarField
        //  onto a new patch
        turbulentIntensityKineticEnergyInletFvPatchScalarField
        (
            const turbulentIntensityKineticEnergyInletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentIntensityKineticEnergyInletFvPatchScalarField
        (
            const turbulentIntensityKineticEnergyInletFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new turbulentIntensityKineticEnergyInletFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        turbulentIntensityKineticEnergyInletFvPatchScalarField
        (
            const turbulentIntensityKineticEnergyInletFvPatchScalarField&,
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
                new turbulentIntensityKineticEnergyInletFvPatchScalarField
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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
