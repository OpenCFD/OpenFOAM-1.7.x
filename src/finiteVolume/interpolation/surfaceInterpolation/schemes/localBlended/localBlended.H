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
    Foam::localBlended

Description
    Two-scheme localBlended differencing scheme.

SourceFiles
    localBlended.C

\*---------------------------------------------------------------------------*/

#ifndef localBlended_H
#define localBlended_H

#include "surfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class localBlended Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class localBlended
:
    public surfaceInterpolationScheme<Type>
{
    // Private Member Functions

        //- Scheme 1
        tmp<surfaceInterpolationScheme<Type> > tScheme1_;

        //- Scheme 2
        tmp<surfaceInterpolationScheme<Type> > tScheme2_;


        //- Disallow default bitwise copy construct
        localBlended(const localBlended&);

        //- Disallow default bitwise assignment
        void operator=(const localBlended&);


public:

    //- Runtime type information
    TypeName("localBlended");


    // Constructors

        //- Construct from mesh and Istream. 
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        localBlended
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, is)
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        localBlended
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            surfaceInterpolationScheme<Type>(mesh),
            tScheme1_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            ),
            tScheme2_
            (
                surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
            )
        {}


    // Member Functions

        //- Return the interpolation weighting factors
        tmp<surfaceScalarField> weights
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const
        {
            const surfaceScalarField& blendingFactor =
                this->mesh().objectRegistry::
                lookupObject<const surfaceScalarField>
                (
                    word(vf.name() + "BlendingFactor")
                );

            return
                blendingFactor*tScheme1_().weights(vf)
              + (scalar(1) - blendingFactor)*tScheme2_().weights(vf);
        }

        //- Return the face-interpolate of the given cell field
        //  with explicit correction
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        interpolate(const GeometricField<Type, fvPatchField, volMesh>& vf) const
        {
            const surfaceScalarField& blendingFactor = 
            (
                this->mesh().objectRegistry::
                lookupObject<const surfaceScalarField> 
                (
                    word(vf.name() + "BlendingFactor")
                )
            );

            return
                blendingFactor*tScheme1_().interpolate(vf)
              + (scalar(1) - blendingFactor)*tScheme2_().interpolate(vf);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
