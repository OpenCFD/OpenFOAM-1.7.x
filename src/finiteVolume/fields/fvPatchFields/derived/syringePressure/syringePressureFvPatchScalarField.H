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
    Foam::syringePressureFvPatchScalarField

Description
    Pressure boundary condition obtained from a 0-D model of the cylinder
    of a syringe.

    The syringe cylinder is defined by its initial volume, piston area and
    velocity profile specified by regions of constant acceleration, speed
    and deceleration.  The gas in the cylinder is described by its initial
    pressure and compressibility which is assumed constant, i.e. isothermal
    expansion/compression.

    Example of the BC specification:
    @verbatim
        outlet
        {
            type            syringePressure;
            Ap              1.388e-6; // syringe piston area
            Sp              0.01;     // syringe piston speed
            VsI             1.388e-8; // Initial syringe volume
            tas             0.001;    // Start of piston acceleration
            tae             0.002;    // End of piston acceleration
            tds             0.005;    // Start of piston deceleration
            tde             0.006;    // end of piston deceleration
            psI             1e5;      // Initial syringe pressure
            psi             1e-5;     // Gas compressibility
            ams             0;        // Added (or removed) gas mass
                                      // Initially 0 but used for restarting.
        }
    @endverbatim

SourceFiles
    syringePressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef syringePressureFvPatchScalarField_H
#define syringePressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class syringePressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class syringePressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- syringe piston area
        scalar Ap_;

        //- syringe piston speed
        scalar Sp_;

        //- Initial syringe volume
        scalar VsI_;

        //- Start of piston acceleration
        scalar tas_;

        //- End of piston acceleration
        scalar tae_;

        //- Start of piston deceleration
        scalar tds_;

        //- end of piston deceleration
        scalar tde_;

        //- Initial syringe pressure
        scalar psI_;

        //- Gas compressibility
        scalar psi_;

        //- Added gas mass
        scalar ams_;

        //- Added gas mass at previous time step
        scalar ams0_;

        //- Current time index used to store ms0_
        label curTimeIndex_;


        //- return the volume of the syringe at time t
        scalar Vs(const scalar t) const;


public:

    //- Runtime type information
    TypeName("syringePressure");


    // Constructors

        //- Construct from patch and internal field
        syringePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        syringePressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given syringePressureFvPatchScalarField
        //  onto a new patch
        syringePressureFvPatchScalarField
        (
            const syringePressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new syringePressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy
        syringePressureFvPatchScalarField
        (
            const syringePressureFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        syringePressureFvPatchScalarField
        (
            const syringePressureFvPatchScalarField&,
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
                new syringePressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

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
