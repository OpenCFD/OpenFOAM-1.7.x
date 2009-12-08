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
    Foam::sixDoFRigidBodyDisplacementPointPatchVectorField

Description
    Foam::sixDoFRigidBodyDisplacementPointPatchVectorField

SourceFiles
    sixDoFRigidBodyDisplacementPointPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef sixDoFRigidBodyDisplacementPointPatchVectorField_H
#define sixDoFRigidBodyDisplacementPointPatchVectorField_H

#include "fixedValuePointPatchField.H"
#include "sixDoFRigidBodyMotion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class sixDoFRigidBodyDisplacementPointPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class sixDoFRigidBodyDisplacementPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
    // Private data

        //- Six dof motion object
        sixDoFRigidBodyMotion motion_;

        //- Reference positions of points on the patch
        pointField p0_;

        //- Reference density required by the forces object for
        //  incompressible calculations
        scalar rhoInf_;


public:

    //- Runtime type information
    TypeName("sixDoFRigidBodyDisplacement");


    // Constructors

        //- Construct from patch and internal field
        sixDoFRigidBodyDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sixDoFRigidBodyDisplacementPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<vector> onto a new patch
        sixDoFRigidBodyDisplacementPointPatchVectorField
        (
            const sixDoFRigidBodyDisplacementPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<vector> > clone() const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new sixDoFRigidBodyDisplacementPointPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        sixDoFRigidBodyDisplacementPointPatchVectorField
        (
            const sixDoFRigidBodyDisplacementPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector> > clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector> >
            (
                new sixDoFRigidBodyDisplacementPointPatchVectorField
                (
                    *this,
                    iF
                )
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
