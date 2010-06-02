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
    Foam::octreeDataPoint

Description
    Encapsulation of data needed for octree searches.

    Used for searching for nearest point. No bounding boxes around points.
    Only overlaps and calcNearest are implemented, rest makes little sense.
    Holds (reference to) pointField.

SourceFiles
    octreeDataPoint.C
    octreeDataPointTreaLeaf.H   (template specialization of treeleaf)
    octreeDataPointTreeLeaf.C   (template specialization of treeleaf)

\*---------------------------------------------------------------------------*/

#ifndef octreeDataPoint_H
#define octreeDataPoint_H

#include "point.H"
#include "pointField.H"
#include "treeBoundBox.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class octree;

/*---------------------------------------------------------------------------*\
                           Class octreeDataPoint Declaration
\*---------------------------------------------------------------------------*/

class octreeDataPoint
{
    // Private data

        const pointField& points_;

public:

    // Constructors

        //- Construct from components. Holds reference to points!
        explicit octreeDataPoint(const pointField&);


    // Member Functions

        // Access

            const pointField& points() const
            {
                return points_;
            }

            label size() const
            {
                return points_.size();
            }

        // Search

            //- Get type of sample
            label getSampleType
            (
                const octree<octreeDataPoint>&,
                const point&
            ) const;

            //- Does (bb of) shape at index overlap bb
            bool overlaps
            (
                const label index,
                const treeBoundBox& sampleBb
            ) const;

            //- Does shape at index contain sample
            bool contains
            (
                const label index,
                const point& sample
            ) const;

            //- Segment (from start to end) intersection with shape
            //  at index. If intersects returns true and sets intersectionPoint
            bool intersects
            (
                const label index,
                const point& start,
                const point& end,
                point& intersectionPoint
            ) const;

            //- Sets newTightest to bounding box (and returns true) if
            //  nearer to sample than tightest bounding box. Otherwise
            //  returns false.
            bool findTightest
            (
                const label index,
                const point& sample,
                treeBoundBox& tightest
            ) const;

            //- Given index get unit normal and calculate (numerical) sign
            //  of sample.
            //  Used to determine accuracy of calcNearest or inside/outside.
            //  Note: always returns GREAT since no inside/outside.
            scalar calcSign
            (
                const label index,
                const point& sample,
                vector& n
            ) const;


            //- Calculates nearest (to sample) point on/in shape.
            //  Returns point and mag(nearest - sample)
            scalar calcNearest
            (
                const label index,
                const point& sample,
                point& nearest
            ) const;

            //- Calculates nearest (to line segment) point in shape.
            //  Returns distance and both point.
            scalar calcNearest
            (
                const label index,
                const linePointRef& ln,
                point& linePt,          // nearest point on line
                point& shapePt          // nearest point on shape
            ) const;



        // Write

            //- Write shape at index
            void write(Ostream& os, const label index) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const octreeDataPoint&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
