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
    Foam::volPointInterpolation

Description
    Foam::volPointInterpolation

SourceFiles
    volPointInterpolation.C
    volPointInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef volPointInterpolation_H
#define volPointInterpolation_H

#include "MeshObject.H"
#include "pointPatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class pointMesh;

/*---------------------------------------------------------------------------*\
                       Class volPointInterpolation Declaration
\*---------------------------------------------------------------------------*/

class volPointInterpolation
:
    public MeshObject<fvMesh, volPointInterpolation>
{
    // Private data

        //- Boundary interpolation engine.
        pointPatchInterpolation boundaryInterpolator_;

        //- Interpolation scheme weighting factor array.
        scalarListList pointWeights_;


    // Private member functions

        //- Construct point weighting factors
        void makeWeights();

        //- Disallow default bitwise copy construct
        volPointInterpolation(const volPointInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const volPointInterpolation&);


public:

    // Declare name of the class and its debug switch
    ClassName("volPointInterpolation");


    // Constructors

        //- Constructor given fvMesh and pointMesh.
        explicit volPointInterpolation(const fvMesh&);


    // Destructor

        ~volPointInterpolation();


    // Member functions

        // Access

            const fvMesh& mesh() const
            {
                return boundaryInterpolator_.mesh();
            }


        // Edit

            //- Update mesh topology using the morph engine
            void updateMesh();

            //- Correct weighting factors for moving mesh.
            bool movePoints();


    // Interpolation functions

        //- Interpolate internal field from volField to pointField
        //  using inverse distance weighting
        template<class Type>
        void interpolateInternalField
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Interpolate from volField to pointField
        //  using inverse distance weighting
        template<class Type>
        void interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Interpolate volField using inverse distance weighting
        //  returning pointField with the same patchField types
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            const wordList& patchFieldTypes
        ) const;

        //- Interpolate tmp<volField> using inverse distance weighting
        //  returning pointField with the same patchField types
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh> >&,
            const wordList& patchFieldTypes
        ) const;

        //- Interpolate volField using inverse distance weighting 
        //  returning pointField
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Interpolate tmp<volField> using inverse distance weighting 
        //  returning pointField
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh> > interpolate
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh> >&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "volPointInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
