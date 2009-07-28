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
    spline class : define a basic spline on nKnots knots - this
    does not go anywhere near these knots (will act as a base type for
    various splines that will have real uses)

\*---------------------------------------------------------------------------*/

#include "spline.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
spline::spline(const pointField& a)
:
    knots_(a)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// function B : this is the blending function for constructing the
// spline

scalar spline::B(const scalar& tau) const
{
    scalar value = 0.0;

    if (tau>=2.0 || tau<=-2.0)
    {
        value = 0.0;
    }
    else if (tau<=-1.0)
    {
        value = pow((2.0 + tau),3.0)/6.0;
    }
    else if (tau<=0.0)
    {
        value = (4.0 - 6.0*tau*tau - 3.0*tau*tau*tau)/6.0;
    }
    else if (tau<=1.0)
    {
        value = (4.0 - 6.0*tau*tau + 3.0*tau*tau*tau)/6.0;
    }
    else if (tau<=2.0)
    {
        value = pow((2.0 - tau),3.0)/6.0;
    }
    else
    {
        FatalErrorIn("spline::B(const scalar&)")
            << "How the hell did we get here???, "
            << "tau = " << tau
            << abort(FatalError);
    }

    return value;
}


// position : returns the position along the spline corresponding to the
// variable mu1

vector spline::position(const scalar mu1) const
{
    vector tmp(vector::zero);

    for (register label i=0; i<knots_.size(); i++)
    {
        tmp += B((knots_.size() - 1)*mu1 - i)*knots_[i];
    }

    return tmp;
}


//- Return the length of the spline curve
scalar spline::length() const
{
    notImplemented("spline::length() const");
    return 1.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
