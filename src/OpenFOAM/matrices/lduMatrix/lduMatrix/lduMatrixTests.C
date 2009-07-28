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

Description
    Convergence and singularity tests for solvers.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::lduMatrix::solverPerformance::checkConvergence
(
    const scalar Tolerance,
    const scalar RelTolerance
)
{
    if (debug >= 2)
    {
        Info<< solverName_
            << ":  Iteration " << noIterations_
            << " residual = " << finalResidual_
            << endl;
    }

    if
    (
        finalResidual_ < Tolerance
     || (
            RelTolerance > SMALL
         && finalResidual_ <= RelTolerance*initialResidual_
        )
  // || (solverName == "symSolve" && iter == 0)
    )
    {
        converged_ = true;
    }
    else
    {
        converged_ = false;
    }

    return converged_;
}


bool Foam::lduMatrix::solverPerformance::checkSingularity
(
    const scalar residual
)
{
    if (residual > VSMALL)
    {
        singular_ = false;
    }
    else
    {
        singular_ = true;
    }

    return singular_;
}


void Foam::lduMatrix::solverPerformance::print() const
{
    if (debug)
    {
        Info<< solverName_ << ":  Solving for " << fieldName_;

        if (singular())
        {
            Info<< ":  solution singularity" << endl;
        }
        else
        {
            Info<< ", Initial residual = " << initialResidual_
                << ", Final residual = " << finalResidual_
                << ", No Iterations " << noIterations_
                << endl;
        }
    }
}


// ************************************************************************* //
