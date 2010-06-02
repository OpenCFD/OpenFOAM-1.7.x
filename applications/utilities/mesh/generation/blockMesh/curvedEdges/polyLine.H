/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::polyLine

Description
    A series of straight line segments, which can also be interpreted as
    a series of control points for splines, etc.

ToDo
    A future implementation could also handle a closed polyLine.

SourceFiles
    polyLine.C

\*---------------------------------------------------------------------------*/

#ifndef polyLine_H
#define polyLine_H

#include "pointField.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class polyLine Declaration
\*---------------------------------------------------------------------------*/


class polyLine
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        polyLine(const polyLine&);

        //- Disallow default bitwise assignment
        void operator=(const polyLine&);

protected:

    // Protected data

        //- The control points or ends of each segments
        pointField points_;

        //- The real line length
        scalar lineLength_;

        //- The rational (0-1) cumulative parameter value for each point
        scalarList param_;

    // Protected member functions

        //- Precalculate the rational cumulative parameter value
        //  and the line-length
        void calcParam();


        //- Return the line segment and the local parameter [0..1]
        //  corresponding to the global lambda [0..1]
        label localParameter(scalar& lambda) const;

public:

    // Constructors

        //- Construct from components
        polyLine
        (
            const pointField&,
            const bool notImplementedClosed = false
        );


    // Member Functions

        //- Return const-access to the control-points
        const pointField& points() const;

        //- Return the number of line segments
        label nSegments() const;

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar) const;

        //- Return the point position corresponding to the local parameter
        //  0 <= lambda <= 1 on the given segment
        point position(const label segment, const scalar) const;

        //- Return the length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
