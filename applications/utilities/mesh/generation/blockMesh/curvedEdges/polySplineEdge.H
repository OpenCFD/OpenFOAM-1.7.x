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

Class
    Foam::polySplineEdge

Description
    polySplineEdge class : representation of a spline via a polyLine

SourceFiles
    polySplineEdge.C

\*---------------------------------------------------------------------------*/

#ifndef polySplineEdge_H
#define polySplineEdge_H

#include "curvedEdge.H"
#include "polyLine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polySplineEdge Declaration
\*---------------------------------------------------------------------------*/

class polySplineEdge
:
    public curvedEdge,
    public polyLine
{
    // Private data

        pointField otherKnots_;


    // Private member functions

        pointField intervening
        (
            const pointField& otherknots,
             const label nbetweenKnots,
            const vector&, const vector&
        );

        label nsize(const label otherknotsSize, const label nbetweenKnots)
        {
            return otherknotsSize*(1 + nbetweenKnots) + nbetweenKnots + 2;
        }


public:

    //- Runtime type information
    TypeName("polySpline");


    // Constructors

        //- Construct from components
        polySplineEdge
        (
            const pointField& ps,
            const label start,
            const label end,
            const pointField& otherknots,
            const label nInterKnots = 20
        );

        //- Construct from Istream setting pointsList
        polySplineEdge(const pointField& points, Istream& is);


    // Destructor

        virtual ~polySplineEdge(){}


    // Member Functions

        //- Return the position of a point on the curve given by
        //  the parameter 0 <= lambda <= 1
        vector position(const scalar mu) const;

        //- Return the length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
