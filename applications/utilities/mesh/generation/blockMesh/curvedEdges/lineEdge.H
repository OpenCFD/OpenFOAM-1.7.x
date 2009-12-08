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
    Foam::lineEdge

Description
    A straight edge between the start point and the end point.

SourceFiles
    lineEdge.C

\*---------------------------------------------------------------------------*/

#ifndef lineEdge_H
#define lineEdge_H

#include "curvedEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class lineEdge Declaration
\*---------------------------------------------------------------------------*/


class lineEdge
:
    public curvedEdge
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        lineEdge(const lineEdge&);

        //- Disallow default bitwise assignment
        void operator=(const lineEdge&);

public:

    //- Runtime type information
    TypeName("line");

    // Constructors

        //- Construct from components
        lineEdge(const pointField&, const label start, const label end);

        //- Construct from Istream with a pointField
        lineEdge(const pointField&, Istream&);


        //- Destructor
        virtual ~lineEdge();


    // Member Functions

        //- Return the point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar) const;

        //- Return the length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
