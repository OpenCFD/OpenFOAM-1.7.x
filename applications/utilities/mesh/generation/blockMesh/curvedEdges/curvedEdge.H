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
    Foam::curvedEdge

Description
    curvedEdges : library functions that will define a curvedEdge in space
    parameterised for 0<lambda<1 from the beginning point to the end point.
    This file contains the abstract base class curvedEdge.

SourceFiles
    curvedEdge.C

\*---------------------------------------------------------------------------*/

#ifndef curvedEdges_H
#define curvedEdges_H

#include "pointField.H"
#include "typeInfo.H"
#include "HashTable.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class curvedEdge Declaration
\*---------------------------------------------------------------------------*/

class curvedEdge
{
protected:

    // Protected data

        const pointField& points_;
        const label start_, end_;

public:

    //- Runtime type information
    TypeName("curvedEdge");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            curvedEdge,
            Istream,
            (
                const pointField& points,
                Istream& is
            ),
            (points, is)
        );


    // Constructors

        //- Construct from components
        curvedEdge
        (
            const pointField& points,
            const label start,
            const label end
        );

        //- Construct from Istream setting pointsList
        curvedEdge(const pointField&, Istream&);

        //- Copy construct
        curvedEdge(const curvedEdge&);

        //- Clone function
        virtual autoPtr<curvedEdge> clone() const;

        //- New function which constructs and returns pointer to a curvedEdge
        static autoPtr<curvedEdge> New(const pointField&, Istream&);


    // Destructor

        virtual ~curvedEdge(){}


    // Member Functions

        //- Return label of start point
        label start() const
        {
            return start_;
        }

        //- Return label of end point
        label end() const
        {
            return end_;
        }

        //- Compare the given start and end points with those of this curve
        bool compare(const label start, const label end) const
        {
            return
            (
                (start_ == start && end_ == end)
             || (start_ == end && end_ == start)
            );
        }

        //- Return the position of a point on the curve given by
        //  the parameter 0 <= lambda <= 1
        virtual vector position(const scalar) const = 0;

        //- Return the length of the curve
        virtual scalar length() const = 0;

    // Member operators

        void operator=(const curvedEdge&);

    // Ostream operator

        friend Ostream& operator<<(Ostream&, const curvedEdge&);
};


//- Return the complete knotList by adding the start and end points to the
//  given list
pointField knotlist
(
    const pointField& points,
    const label start,
    const label end,
    const pointField& otherknots
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
