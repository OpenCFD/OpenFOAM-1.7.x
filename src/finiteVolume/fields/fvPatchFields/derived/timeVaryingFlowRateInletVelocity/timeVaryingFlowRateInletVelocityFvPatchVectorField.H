/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-2010 OpenCFD Ltd.
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
    Foam::timeVaryingFlowRateInletVelocityFvPatchVectorField

Description
    A time-varying form of a flow normal vector boundary condition.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            timeVaryingFlowRateInletVelocity;
        flowRate        0.2;         // Volumetric/mass flow rate [m3/s or kg/s]
        value           uniform (0 0 0);  // placeholder
        fileName        "$FOAM_CASE/time-series";
        outOfBounds     repeat;           // (error|warn|clamp|repeat)
    }
    @endverbatim

Note
    - The value is positive inwards
    - may not work correctly for transonic inlets!
    - strange behaviour with potentialFoam since the U equation is not solved

See Also
    Foam::interpolationTable and Foam::flowRateInletVelocityFvPatchVectorField

SourceFiles
    timeVaryingFlowRateInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingFlowRateInletVelocityFvPatchVectorField_H
#define timeVaryingFlowRateInletVelocityFvPatchVectorField_H

#include "flowRateInletVelocityFvPatchVectorField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
        Class timeVaryingFlowRateInletVelocityFvPatch Declaration
\*---------------------------------------------------------------------------*/

class timeVaryingFlowRateInletVelocityFvPatchVectorField
:
    public flowRateInletVelocityFvPatchVectorField
{
    // Private data

        //- the time series being used, including the bounding treatment
        interpolationTable<scalar> timeSeries_;


public:

    //- Runtime type information
    TypeName("timeVaryingFlowRateInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingFlowRateInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        timeVaryingFlowRateInletVelocityFvPatchVectorField
        (
            const timeVaryingFlowRateInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingFlowRateInletVelocityFvPatchVectorField
        (
            const timeVaryingFlowRateInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new timeVaryingFlowRateInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingFlowRateInletVelocityFvPatchVectorField
        (
            const timeVaryingFlowRateInletVelocityFvPatchVectorField&,
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
                new timeVaryingFlowRateInletVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        // Access

            //- Return the time series used
            const interpolationTable<scalar>& timeSeries() const
            {
                return timeSeries_;
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
