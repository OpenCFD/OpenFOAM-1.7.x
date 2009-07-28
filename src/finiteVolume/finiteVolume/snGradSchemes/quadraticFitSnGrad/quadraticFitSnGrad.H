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
    quadraticFitSnGrad

Description
    Simple central-difference snGrad scheme with quadratic fit correction from
    a larger stencil.

SourceFiles
    quadraticFitSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticFitSnGrad_H
#define quadraticFitSnGrad_H

#include "snGradScheme.H"
#include "quadraticFitSnGradData.H"
#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class quadraticFitSnGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class quadraticFitSnGrad
:
    public snGradScheme<Type>
{
    // Private Data
        //- weights for central stencil
        const scalar centralWeight_;

    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const quadraticFitSnGrad&);


public:

    //- Runtime type information
    TypeName("quadraticFit");


    // Constructors

        //- Construct from mesh and scheme data
        quadraticFitSnGrad
        (
            const fvMesh& mesh,
            const scalar centralWeight
        )
        :
            snGradScheme<Type>(mesh),
            centralWeight_(centralWeight)
        {}


        //- Construct from mesh and data stream
        quadraticFitSnGrad(const fvMesh& mesh, Istream& is)
        :
            snGradScheme<Type>(mesh),
            centralWeight_(readScalar(is))
        {}


    // Destructor

        virtual ~quadraticFitSnGrad() {}


    // Member Functions

        //- Return the interpolation weighting factors for the given field
        virtual tmp<surfaceScalarField> deltaCoeffs
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const
        {
            return this->mesh().deltaCoeffs();
        }

        //- Return true if this scheme uses an explicit correction
        virtual bool corrected() const
        {
            return true;
        }

        //- Return the explicit correction to the quadraticFitSnGrad
        //  for the given field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
        correction(const GeometricField<Type, fvPatchField, volMesh>& vf) const
        {
            const fvMesh& mesh = this->mesh();

            const quadraticFitSnGradData& qfd = quadraticFitSnGradData::New
            (
                mesh,
                centralWeight_
            );

            const extendedCellToFaceStencil& stencil = qfd.stencil();
            const List<scalarList>& f = qfd.fit();

            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > sft
                = stencil.weightedSum(vf, f);

            sft().dimensions() /= dimLength;

            return sft;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
