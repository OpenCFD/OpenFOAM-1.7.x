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
    Foam::octreeDataTriSurface

Description
    Encapsulates data for octree searches on triSurface.

SourceFiles
    octreeDataTriSurface.C

\*---------------------------------------------------------------------------*/

#ifndef octreeDataTriSurface_H
#define octreeDataTriSurface_H

#include "treeBoundBoxList.H"
#include "labelList.H"
#include "point.H"
#include "triSurface.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class Type> class octree;

/*---------------------------------------------------------------------------*\
                           Class octreeDataTriSurface Declaration
\*---------------------------------------------------------------------------*/

class octreeDataTriSurface
{
    // Static data

        //- tolerance on linear dimensions
        static scalar tol;

    // Private data

        const triSurface& surface_;

        const treeBoundBoxList allBb_;

        // Extra data to speed up distance searches.
        // Triangles expressed as base + spanning vectors
        pointField base_;
        pointField E0_;
        pointField E1_;
        scalarList a_;
        scalarList b_;
        scalarList c_;


    // Private Static Functions

        //- fast triangle nearest point calculation. Returns point in E0, E1
        //  coordinate system:  base + s*E0 + t*E1
        static void nearestCoords
        (
            const point& base,
            const point& E0,
            const point& E1,
            const scalar a,
            const scalar b,
            const scalar c,
            const point& P,
            scalar& s,
            scalar& t
        );

        //- Calculate bounding boxes for triangles
        static treeBoundBoxList calcBb(const triSurface&);

    // Private Member Functions

        //- nearest point in xyz coord system
        point nearestPoint(const label index, const point& P) const;

public:

    // Declare name of the class and its debug switch
    ClassName("octreeDataTriSurface");


    // Constructors

        //- Construct from triSurface. Holds reference. Bounding box
        //  calculated from triangle points.
        octreeDataTriSurface(const triSurface&);

        //- Construct from triSurface and bounding box.
        //  Holds references.
        octreeDataTriSurface(const triSurface&, const treeBoundBoxList&);


    // Member Functions

        // Access

            const triSurface& surface() const
            {
                return surface_;
            }

            const treeBoundBoxList& allBb() const
            {
                return allBb_;
            }

            label size() const
            {
                return allBb_.size();
            }

        // Search

            //- Get type of sample
            label getSampleType
            (
                const octree<octreeDataTriSurface>&,
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
            scalar calcSign
            (
                const label index,
                const point& sample,
                vector& n
            ) const;

            //- Calculates nearest (to sample) point in shape.
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

            // Write shape at index
            void write(Ostream& os, const label index) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, octreeDataTriSurface&);
        friend Ostream& operator<<(Ostream&, const octreeDataTriSurface&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
