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
    BSpline : cubic spline going through all the knots

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "BSpline.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

pointField BSpline::findKnots
(
    const pointField& allknots,
    const vector& fstend,
    const vector& sndend
)
{
    label newnKnots(allknots.size() + 2);
    label NKnots(allknots.size());
    pointField newknots(newnKnots);

    // set up 1/6 and 2/3 which are the matrix elements throughout most
    // of the matrix

    register scalar oneSixth = 1.0/6.0;
    register scalar twoThird = 2.0/3.0;

    simpleMatrix<vector> M(newnKnots);

    // set up the matrix

    M[0][0] = -0.5*scalar(NKnots - 1);
    M[0][2] =  0.5*scalar(NKnots - 1);

    for (register label i=1; i<newnKnots-1; i++)
    {
        M[i][i-1] = oneSixth;
        M[i][i] = twoThird;
        M[i][i+1] = oneSixth;
    }

    M[newnKnots - 1][newnKnots - 3] = -0.5*scalar(NKnots - 1);
    M[newnKnots - 1][newnKnots - 1] =  0.5*scalar(NKnots - 1);

    // set up the vector

    for (label i=1; i<=NKnots; i++)
    {
        M.source()[i] = allknots[i-1];
    }

    // set the gradients at the two ends

    if (mag(fstend)<1e-8)
    {
        // set to the default : forward differences on the end knots
        M.source()[0] = allknots[1] - allknots[0];
        M.source()[0] /= mag(M.source()[0]);

        M.source()[NKnots+1] = M.source()[NKnots-1] - M.source()[NKnots];
        M.source()[NKnots+1] /= mag(M.source()[NKnots+1]);
    }
    else
    {
        // set to the gradient vectors provided
        M.source()[0] = fstend/mag(fstend);
        M.source()[NKnots+1] = sndend/mag(sndend);
    }

    // invert the equation to find the control knots

    newknots = M.solve();

    return newknots;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
BSpline::BSpline(const pointField& Knots)
:
    spline(findKnots(Knots))
{}


// Construct from components
BSpline::BSpline
(
    const pointField& Knots,
    const vector& fstend,
    const vector& sndend
)
:
    spline(findKnots(Knots, fstend, sndend))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the real position of a point on the curve given by
//  the parameter 0 <= lambda <= 1
vector BSpline::realPosition(scalar mu)
{
    return spline::position(mu);
}


//- Return the position of a point on the curve given by
//  the parameter 0 <= lambda <= 1
vector BSpline::position(const scalar mu) const
{
    return spline::position((1.0/(nKnots() - 1))*(1.0 + mu*(nKnots() - 3)));
}


//- Return the length of the curve
scalar BSpline::length() const
{
    notImplemented("BSpline::length() const");
    return 1.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
