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
    Foam::fvMotionSolver

Description
    Virtual base class for finite volume mesh motion solvers.

    The boundary motion is set as a boundary condition on the motion velocity
    variable motionU.

SourceFiles
    fvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef fvMotionSolver_H
#define fvMotionSolver_H

#include "motionSolver.H"
#include "pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class fvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class fvMotionSolver
:
    public motionSolver
{

protected:

    // Protected data

        //- The fvMesh to be moved
        const fvMesh& fvMesh_;


    // Protected member functions

        //- Create the corresponding patch types for cellMotion from those
        //  of the given pointMotion
        template<class Type>
        wordList cellMotionBoundaryTypes
        (
            const typename GeometricField<Type, pointPatchField, pointMesh>::
            GeometricBoundaryField& pmUbf
        ) const;


public:

    //- Runtime type information
    TypeName("fvMotionSolver");


    // Constructors

        //- Construct from polyMesh
        fvMotionSolver(const polyMesh& mesh);


    // Destructor

        virtual ~fvMotionSolver();


    // Member Functions

        //- Return reference to the fvMesh to be moved
        const fvMesh& mesh() const
        {
            return fvMesh_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const = 0;

        //- Solve for motion
        virtual void solve() = 0;

        //- Update the pointMesh corresponding to the new points
        void movePoints(const pointField&);

        //- Update the mesh corresponding to given map
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fvMotionSolverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
