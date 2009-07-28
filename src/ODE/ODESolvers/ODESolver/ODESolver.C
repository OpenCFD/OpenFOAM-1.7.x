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

\*---------------------------------------------------------------------------*/

#include "ODESolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::ODESolver, 0);
namespace Foam
{
    defineRunTimeSelectionTable(ODESolver, ODE);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ODESolver::ODESolver(const ODE& ode)
:
    n_(ode.nEqns()),
    yScale_(n_),
    dydx_(n_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ODESolver::solve
(
    const ODE& ode,
    const scalar xStart,
    const scalar xEnd,
    scalarField& y,
    const scalar eps,
    scalar& hEst
) const
{
    const label MAXSTP = 10000;

    scalar x = xStart;
    scalar h = hEst;

    for (label nStep=0; nStep<MAXSTP; nStep++)
    {
        ode.derivatives(x, y, dydx_);

        for (label i=0; i<n_; i++)
        {
            yScale_[i] = mag(y[i]) + mag(dydx_[i]*h) + SMALL;
        }

        if ((x + h - xEnd)*(x + h - xStart) > 0.0)
        {
            h = xEnd - x;
        }

        scalar hNext, hDid;
        solve(ode, x, y, dydx_, eps, yScale_, h, hDid, hNext);

        if ((x - xEnd)*(xEnd - xStart) >= 0.0)
        {
            hEst = hNext;
            return;
        }

        h = hNext;
    }

    FatalErrorIn
    (
        "ODESolver::solve"
        "(const ODE& ode, const scalar xStart, const scalar xEnd,"
        "scalarField& yStart, const scalar eps, scalar& hEst) const"
    )   << "Too many integration steps"
        << exit(FatalError);
}

// ************************************************************************* //
