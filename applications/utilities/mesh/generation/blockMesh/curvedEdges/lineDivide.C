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
    lineDivide class : divides a line into segments

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "lineDivide.H"
#include "curvedEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
lineDivide::lineDivide(const curvedEdge& bc, const label n, const scalar xratio)
:
    points_(n + 1),
    divisions_(n + 1),
    noPoints_(n)
{
    scalar np(n);
    scalar lambda(0.0);

    if (xratio == 1.0)
    {
        scalar y(1.0/np);
        for (label i=0; i<=noPoints_; i++)
        {
            lambda = scalar(i)/np;
            points_[i] = bc.position(lambda);
            divisions_[i] = y*i;
        }
    }
    else
    {
        points_[0] = bc.position(0.0);
        divisions_[0] = 0.0;
        scalar xrpower = 1.0;

        for (label i=1; i<=noPoints_; i++)
        {
            lambda = (1.0 - pow(xratio, i))/(1.0 - pow(xratio, np));
            points_[i] = bc.position(lambda);
            divisions_[i] = lambda;
            xrpower *= xratio;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const pointField& lineDivide::points() const
{
    return points_;
}


const scalarList& lineDivide::lambdaDivisions() const
{
    return divisions_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
