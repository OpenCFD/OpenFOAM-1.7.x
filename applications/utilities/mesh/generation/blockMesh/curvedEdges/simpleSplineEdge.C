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
    simpleSplineEdge : the actual access class for Bspline

\*---------------------------------------------------------------------------*/

#include "simpleSplineEdge.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(simpleSplineEdge, 0);
addToRunTimeSelectionTable(curvedEdge, simpleSplineEdge, Istream);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
simpleSplineEdge::simpleSplineEdge
(
    const pointField& points,
    const label start,
    const label end,
    const pointField& otherknots
)
:
    curvedEdge(points, start, end),
    BSpline(knotlist(points, start, end, otherknots))
{}


// Construct from components
simpleSplineEdge::simpleSplineEdge
(
    const pointField& points,
    const label start,
    const label end,
    const pointField& otherknots,
    const vector& fstend,
    const vector& sndend
)
:
    curvedEdge(points, start, end),
    BSpline(knotlist(points, start, end, otherknots), fstend, sndend)
{}


// Construct from Istream
simpleSplineEdge::simpleSplineEdge(const pointField& points, Istream& is)
:
    curvedEdge(points, is),
    BSpline(knotlist(points, start_, end_, pointField(is)))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the position of a point on the simple spline curve given by
//  the parameter 0 <= lambda <= 1
vector simpleSplineEdge::position(const scalar mu) const
{
    return BSpline::position(mu);
}


//- Return the length of the simple spline curve
scalar simpleSplineEdge::length() const
{
    notImplemented("simpleSplineEdge::length() const");
    return 1.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
