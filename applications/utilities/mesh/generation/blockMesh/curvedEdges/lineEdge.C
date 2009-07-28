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
    line class : defines a straight line between the start point and the
    end point

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "lineEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(lineEdge, 0);

// Add the curvedEdge constructor functions to the hash tables
curvedEdge::addIstreamConstructorToTable<lineEdge>
    addLineEdgeIstreamConstructorToTable_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
lineEdge::lineEdge
(
    const pointField& points,
    const label start,
    const label end
)
:
    curvedEdge(points, start, end),
    startPoint_(points_[start_]),
    direction_(points_[end_] - points_[start_])
{}


// Construct from Istream
lineEdge::lineEdge(const pointField& points, Istream& is)
:
    curvedEdge(points, is),
    startPoint_(points_[start_]),
    direction_(points_[end_] - points_[start_])
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

vector lineEdge::position(const scalar lambda) const
{
    if (lambda < 0 || lambda > 1)
    {
        FatalErrorIn("lineEdge::position(const scalar)")
            << "Parameter out of range, lambda = " << lambda
            << abort(FatalError);
    }

    return startPoint_ + lambda*direction_;
}


//- Return the length of the curve
scalar lineEdge::length() const
{
    return mag(direction_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
