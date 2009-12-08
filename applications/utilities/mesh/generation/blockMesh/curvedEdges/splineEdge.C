/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "splineEdge.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(splineEdge, 0);

    addToRunTimeSelectionTable
    (
        curvedEdge,
        splineEdge,
        Istream
    );

    // compatibility with old names
    addNamedToRunTimeSelectionTable
    (
        curvedEdge,
        splineEdge,
        Istream,
        simpleSpline
    );

    // compatibility with old names
    addNamedToRunTimeSelectionTable
    (
        curvedEdge,
        splineEdge,
        Istream,
        polySpline
    );

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::splineEdge::splineEdge
(
    const pointField& points,
    const label start,
    const label end,
    const pointField& internalPoints
)
:
    curvedEdge(points, start, end),
    CatmullRomSpline(appendEndPoints(points, start, end, internalPoints))
{}


Foam::splineEdge::splineEdge(const pointField& points, Istream& is)
:
    curvedEdge(points, is),
    CatmullRomSpline(appendEndPoints(points, start_, end_, pointField(is)))
{
    token t(is);
    is.putBack(t);

    // compatibility - might also have start/end tangents - discard them
    if (t == token::BEGIN_LIST)
    {
        vector tangent0Ignored(is);
        vector tangent1Ignored(is);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::splineEdge::~splineEdge()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::splineEdge::position(const scalar mu) const
{
    return CatmullRomSpline::position(mu);
}


Foam::scalar Foam::splineEdge::length() const
{
    return CatmullRomSpline::length();
}


// ************************************************************************* //
