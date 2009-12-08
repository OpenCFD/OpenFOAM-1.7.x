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

Class
    Foam::BSpline

Description
    An implementation of B-splines.

    In this implementation, the end tangents are created automatically
    by reflection.

    In matrix form, the @e local interpolation on the interval t=[0..1] is
    described as follows:
    @verbatim
    P(t) = 1/6 * [ t^3 t^2 t 1 ] * [ -1  3 -3  1 ] * [ P-1 ]
                                   [  3 -6  3  0 ]   [ P0 ]
                                   [ -3  0  3  0 ]   [ P1 ]
                                   [  1  4  1  0 ]   [ P2 ]
    @endverbatim

    Where P-1 and P2 represent the neighbouring points or the extrapolated
    end points. Simple reflection is used to automatically create the end
    points.

    The spline is discretized based on the chord length of the individual
    segments. In rare cases (sections with very high curvatures), the
    resulting distribution may be sub-optimal.

SeeAlso
    CatmullRomSpline

ToDo
    A future implementation could also handle closed splines.

SourceFiles
    BSpline.C

\*---------------------------------------------------------------------------*/

#ifndef BSpline_H
#define BSpline_H

#include "polyLine.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BSpline Declaration
\*---------------------------------------------------------------------------*/

class BSpline
:
    public polyLine
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        BSpline(const BSpline&);

        //- Disallow default bitwise assignment
        void operator=(const BSpline&);


public:

    // Constructors

        //- Construct from components
        BSpline
        (
            const pointField& knots,
            const bool notImplementedClosed = false
        );


    // Member Functions

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar lambda) const;

        //- Return the point position corresponding to the local parameter
        //  0 <= lambda <= 1 on the given segment
        point position(const label segment, const scalar lambda) const;

        //- Return the length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
