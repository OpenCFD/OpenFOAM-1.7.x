/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::motionDiffusivity

Description
    Abstract base class for cell-centre mesh motion diffusivity.

SourceFiles
    motionDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef motionDiffusivity_H
#define motionDiffusivity_H

#include "fvMotionSolver.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class motionDiffusivity Declaration
\*---------------------------------------------------------------------------*/

class motionDiffusivity
{
    // Private data

        //- Motion solver reference
        const fvMotionSolver& mSolver_;

public:

    //- Runtime type information
    TypeName("motionDiffusivity");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            motionDiffusivity,
            Istream,
            (
                const fvMotionSolver& mSolver,
                Istream& mdData
            ),
            (mSolver, mdData)
        );


    // Selectors

        //- Select null constructed
        static autoPtr<motionDiffusivity> New
        (
            const fvMotionSolver& mSolver,
            Istream& mdData
        );


    // Constructors

        //- Construct for the given fvMotionSolver
        motionDiffusivity(const fvMotionSolver& mSolver);


    // Destructor

        virtual ~motionDiffusivity();


    // Member Functions

        //- Return reference to the motion solver
        const fvMotionSolver& mSolver() const
        {
            return mSolver_;
        }

        //- Return diffusivity field
        virtual tmp<surfaceScalarField> operator()() const = 0;

        //- Correct the motion diffusivity
        virtual void correct() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
