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
    Foam::uniformDiffusivity

Description
    Uniform uniform finite volume mesh motion diffusivity.

SourceFiles
    uniformDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef uniformDiffusivity_H
#define uniformDiffusivity_H

#include "motionDiffusivity.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uniformDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class uniformDiffusivity
:
    public motionDiffusivity
{

protected:

    // Protected data

        surfaceScalarField faceDiffusivity_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        uniformDiffusivity(const uniformDiffusivity&);

        //- Disallow default bitwise assignment
        void operator=(const uniformDiffusivity&);


public:

    //- Runtime type information
    TypeName("uniform");


    // Constructors

        //- Construct for the given fvMotionSolver and data Istream
        uniformDiffusivity
        (
            const fvMotionSolver& mSolver,
            Istream& mdData
        );


    // Destructor

        virtual ~uniformDiffusivity();


    // Member Functions

        //- Return diffusivity field
        virtual tmp<surfaceScalarField> operator()() const
        {
            return faceDiffusivity_;
        }

        //- Do not correct the motion diffusivity
        virtual void correct()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
