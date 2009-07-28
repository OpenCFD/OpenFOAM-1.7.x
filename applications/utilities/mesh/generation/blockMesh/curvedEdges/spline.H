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
    Foam::spline

Description
    spline class : define a basic spline on nKnots knots - this
    does not go anywhere near these knots (will act as a base type for
    various splines that will have real uses)

SourceFiles
    spline.C

\*---------------------------------------------------------------------------*/

#ifndef spline_H
#define spline_H

#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class spline Declaration
\*---------------------------------------------------------------------------*/

class spline
{
    // Private data

        //- The knots defining the spline
        pointField knots_;


    // Private member functions

        //- Blending function for constructing spline
        scalar B(const scalar&) const;


public:

    // Constructors

        //- Construct from components
        spline(const pointField&);


    // Member Functions

        // Access

            //- Return the number of knots in the spline
            label nKnots() const
            {
                return knots_.size();
            }


        //- Return the position of a point on the curve given by
        //  the parameter 0 <= lambda <= 1
        vector position(const scalar) const;

        //- Return the length of the spline curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
