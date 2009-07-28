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
    arcEdge class : defines the arcEdge of a circle in terms of 3 points on its
    circumference

\*---------------------------------------------------------------------------*/

#include "arcEdge.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(arcEdge, 0);

    // Add the curvedEdge constructor functions to the hash tables
    addToRunTimeSelectionTable(curvedEdge, arcEdge, Istream);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::cylindricalCS Foam::arcEdge::calcAngle()
{
    vector a = p2_ - p1_;
    vector b = p3_ - p1_;

    // find centre of arcEdge
    scalar asqr = a & a;
    scalar bsqr = b & b;
    scalar adotb = a & b;

    scalar denom = asqr*bsqr - adotb*adotb;

    if (mag(denom) < VSMALL)
    {
        FatalErrorIn("cylindricalCS arcEdge::calcAngle()")
            << "Invalid arc definition - are the points co-linear?  Denom ="
            << denom
            << abort(FatalError);
    }

    scalar fact = 0.5*(bsqr - adotb)/denom;

    vector centre = 0.5*a + fact*((a ^ b) ^ a);

    centre += p1_;

    // find position vectors w.r.t. the arcEdge centre
    vector r1(p1_ - centre);
    vector r2(p2_ - centre);
    vector r3(p3_ - centre);

    // find angles
    scalar tmp = (r3&r1)/(mag(r3)*mag(r1));
    angle_ = acos(tmp)*180.0/mathematicalConstant::pi;

    // check if the vectors define an exterior or an interior arcEdge
    if (((r1  ^ r2)&(r1 ^ r3)) < 0.0) angle_ = 360 - angle_;

    vector tempAxis(0.0,0.0,0.0);

    if (angle_ <= 180.0)
    {
        tempAxis = r1 ^ r3;

        if (mag(tempAxis)/(mag(r1)*mag(r3)) < 0.001) tempAxis = r1 ^ r2;
    }
    else
    {
        tempAxis = r3 ^ r1;
    }

    radius_ = mag(r3);

    // set up and return the local coordinate system
    return cylindricalCS("tmpCS", centre, tempAxis, r1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::arcEdge::arcEdge
(
    const pointField& points,
    const label start,
    const label end,
    const vector& P2
)
:
    curvedEdge(points, start, end),
    p1_(points_[start_]),
    p2_(P2),
    p3_(points_[end_]),
    cs_(calcAngle())
{}


// Construct from Istream
Foam::arcEdge::arcEdge(const pointField& points, Istream& is)
:
    curvedEdge(points, is),
    p1_(points_[start_]),
    p2_(is),
    p3_(points_[end_]),
    cs_(calcAngle())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::vector Foam::arcEdge::position(const scalar lambda) const
{
    if (lambda < 0 || lambda > 1)
    {
        FatalErrorIn("arcEdge::position(const scalar lambda) const")
            << "Parameter out of range, lambda = " << lambda
            << abort(FatalError);
    }

    if (lambda < SMALL)
    {
        return p1_;
    }
    else if (lambda > 1-SMALL)
    {
        return p3_;
    }
    else
    {
        return cs_.globalPosition(vector(radius_, lambda*angle_, 0.0));
    }
}


//- Return the length of the curve
Foam::scalar Foam::arcEdge::length() const
{
    return angle_*radius_*mathematicalConstant::pi/180.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
