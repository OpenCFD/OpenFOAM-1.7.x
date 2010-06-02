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
    Foam::pointPatchInterpolation

Description
    Foam::pointPatchInterpolation

SourceFiles
    pointPatchInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef pointPatchInterpolation_H
#define pointPatchInterpolation_H

#include "primitivePatchInterpolation.H"
#include "PtrList.H"
#include "volFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "scalarList.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class pointMesh;

/*---------------------------------------------------------------------------*\
                       Class pointPatchInterpolation Declaration
\*---------------------------------------------------------------------------*/

class pointPatchInterpolation
{
    // Private data

        const fvMesh& fvMesh_;

        //- Primitive patch interpolators
        PtrList<primitivePatchInterpolation> patchInterpolators_;

        //- List of patch-patch edge points that require special treatement
        labelList patchPatchPoints_;

        //- Weights for patch-patch boundary points
        scalarListList patchPatchPointWeights_;

        labelList patchPatchPointConstraintPoints_;
        tensorField patchPatchPointConstraintTensors_;


    // Private member functions

        //- Construct addressing for patch-patch boundary points
        void makePatchPatchAddressing();

        //- Construct weights for patch-patch boundary points
        void makePatchPatchWeights();


        //- Disallow default bitwise copy construct
        pointPatchInterpolation(const pointPatchInterpolation&);

        //- Disallow default bitwise assignment
        void operator=(const pointPatchInterpolation&);


public:

    // Declare name of the class and its debug switch
    ClassName("pointPatchInterpolation");


    // Constructors

        //- Constructor given fvMesh and pointMesh.
        pointPatchInterpolation(const fvMesh&);


    // Destructor

        ~pointPatchInterpolation();


    // Member functions

        // Access

            const fvMesh& mesh() const
            {
                return fvMesh_;
            }


        // Edit

            //- Update mesh topology using the morph engine
            void updateMesh();

            //- Correct weighting factors for moving mesh.
            bool movePoints();


    // Interpolation functions

        template<class Type>
        void interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            GeometricField<Type, pointPatchField, pointMesh>&,
            bool overrideFixedValue
        ) const;

        template<class Type>
        void applyCornerConstraints
        (
            GeometricField<Type, pointPatchField, pointMesh>& pf
        ) const;
};


template<>
void pointPatchInterpolation::applyCornerConstraints<scalar>
(
    GeometricField<scalar, pointPatchField, pointMesh>& pf
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "pointPatchInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
