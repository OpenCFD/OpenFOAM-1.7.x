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
    Foam::timeVaryingUniformInletOutletFvPatchField

Description
    A time-varying form of a uniform fixed value boundary condition.

    Example of the boundary condition specification:
    @verbatim
    inlet
    {
        type            timeVaryingUniformInletOutlet;
        fileName        "$FOAM_CASE/time-series";
        outOfBounds     clamp;           // (error|warn|clamp|repeat)
    }
    @endverbatim

Note
    This class is derived directly from a inletOutlet patch rather than from
    a uniformInletOutlet patch.

See Also
    Foam::interpolationTable and Foam::inletOutletFvPatchField

SourceFiles
    timeVaryingUniformInletOutletFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef timeVaryingUniformInletOutletFvPatchField_H
#define timeVaryingUniformInletOutletFvPatchField_H

#include "inletOutletFvPatchField.H"
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class timeVaryingUniformInletOutletFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class timeVaryingUniformInletOutletFvPatchField
:
    public inletOutletFvPatchField<Type>
{
    // Private data

        //- The time series being used, including the bounding treatment
        interpolationTable<Type> timeSeries_;


public:

    //- Runtime type information
    TypeName("timeVaryingUniformInletOutlet");


    // Constructors

        //- Construct from patch and internal field
        timeVaryingUniformInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        timeVaryingUniformInletOutletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        timeVaryingUniformInletOutletFvPatchField
        (
            const timeVaryingUniformInletOutletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        timeVaryingUniformInletOutletFvPatchField
        (
            const timeVaryingUniformInletOutletFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new timeVaryingUniformInletOutletFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        timeVaryingUniformInletOutletFvPatchField
        (
            const timeVaryingUniformInletOutletFvPatchField<Type>&,
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
                new timeVaryingUniformInletOutletFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the time series used
            const interpolationTable<Type>& timeSeries() const
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

#ifdef NoRepository
#   include "timeVaryingUniformInletOutletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
