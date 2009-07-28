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
    Foam::arcEdge

Description
    arcEdge class : defines the arcEdge of a circle in terms of 3 points on its
    circumference

SourceFiles
    arcEdge.C

\*---------------------------------------------------------------------------*/

#ifndef arcEdge_H
#define arcEdge_H

#include "curvedEdge.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class arcEdge Declaration
\*---------------------------------------------------------------------------*/

class arcEdge
:
    public curvedEdge
{
    // Private data

        vector p1_, p2_, p3_;
        scalar angle_;
        scalar radius_;
        cylindricalCS cs_;

        cylindricalCS calcAngle();

public:

    //- Runtime type information
    TypeName("arc");


    // Constructors

        //- Construct from components
        arcEdge
        (
            const pointField& points,
            const label start, const label end,
            const vector&
        );

        //- Construct from Istream setting pointsList
        arcEdge(const pointField& points, Istream&);


    // Destructor

        virtual ~arcEdge(){}


    // Member Functions

        //- Return the position of a point on the curve given by
        //  the parameter 0 <= lambda <= 1
        vector position(const scalar) const;

        //- Return the length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
