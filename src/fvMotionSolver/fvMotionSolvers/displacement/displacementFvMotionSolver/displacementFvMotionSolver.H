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
    Foam::displacementFvMotionSolver.H

Description
    Base class for fvMotionSolvers which calculate displacement.

SourceFiles
    displacementFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementFvMotionSolver_H
#define displacementFvMotionSolver_H

#include "fvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class displacementFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementFvMotionSolver
:
    public fvMotionSolver
{
    // Private data

        //- Reference point field
        pointField points0_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        displacementFvMotionSolver
        (
            const displacementFvMotionSolver&
        );

        //- Disallow default bitwise assignment
        void operator=(const displacementFvMotionSolver&);


public:

    //- Runtime type information
    TypeName("displacementInterpolation");


    // Constructors

        //- Construct from polyMesh and data stream
        displacementFvMotionSolver
        (
            const polyMesh&,
            Istream& msDataUnused
        );


    // Destructor

        ~displacementFvMotionSolver();


    // Member Functions

        //- Return reference to the reference field
        const pointField& points0() const
        {
            return points0_;
        }

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
