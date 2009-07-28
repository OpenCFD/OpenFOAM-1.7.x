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
    Foam::surfaceInterpolation

Description
    Cell to surface interpolation scheme. Included in fvMesh.

SourceFiles
    surfaceInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceInterpolation_H
#define surfaceInterpolation_H

#include "tmp.H"
#include "scalar.H"
#include "fvSchemes.H"
#include "fvSolution.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class surfaceInterpolation Declaration
\*---------------------------------------------------------------------------*/

class surfaceInterpolation
:
    public fvSchemes,
    public fvSolution
{
    // Private data

        // Reference to fvMesh
        const fvMesh& mesh_;

        // Demand-driven data

            //- Central-differencing weighting factors
            mutable surfaceScalarField* weightingFactors_;

            //- Face-gradient difference factors
            mutable surfaceScalarField* differenceFactors_;

            //- Is mesh orthogonal
            mutable bool orthogonal_;

            //- Non-orthogonality correction vectors
            mutable surfaceVectorField* correctionVectors_;


    // Private member functions

        //- Construct central-differencing weighting factors
        void makeWeights() const;

        //- Construct face-gradient difference factors
        void makeDeltaCoeffs() const;

        //- Construct non-orthogonality correction vectors
        void makeCorrectionVectors() const;


protected:

    // Protected member functions

        // Storage management

            //- Clear all geometry and addressing
            void clearOut();


public:

    // Declare name of the class and its debug switch
    ClassName("surfaceInterpolation");


    // Constructors

        //- Construct given an fvMesh
        explicit surfaceInterpolation(const fvMesh&);


    // Destructor

        ~surfaceInterpolation();


    // Member functions

        //- Return reference to weighting factors array
        const surfaceScalarField& weights() const;

        //- Return reference to difference factors array
        const surfaceScalarField& deltaCoeffs() const;

        //- Return whether mesh is orthogonal or not
        bool orthogonal() const;

        //- Return reference to non-orthogonality correction vectors array
        const surfaceVectorField& correctionVectors() const;

        //- Do what is neccessary if the mesh has moved
        bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
