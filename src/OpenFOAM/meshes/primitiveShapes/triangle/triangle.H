/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::triangle

Description
    A triangle primitive used to calculate face normals and swept volumes.

SourceFiles
    triangleI.H

\*---------------------------------------------------------------------------*/

#ifndef triangle_H
#define triangle_H

#include "intersection.H"
#include "vector.H"
#include "pointHit.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// Forward declaration of friend functions and operators

template<class Point, class PointRef> class triangle;

template<class Point, class PointRef>
inline Istream& operator>>
(
    Istream&,
    triangle<Point, PointRef>&
);

template<class Point, class PointRef>
inline Ostream& operator<<
(
    Ostream&,
    const triangle<Point, PointRef>&
);


/*---------------------------------------------------------------------------*\
                           class triangle Declaration
\*---------------------------------------------------------------------------*/

template<class Point, class PointRef>
class triangle
{
    // Private data

        PointRef a_, b_, c_;

    // Private member functions

        //- Fast distance to triangle calculation. From
        //  "Distance Between Point and Trangle in 3D"
        //  David Eberly, Magic Software Inc. Aug. 2002.
        //  Works on function Q giving distance to point and tries to
        //  minimize this.
        static pointHit nearestPoint
        (
            const Point& baseVertex,
            const vector& E0,
            const vector& E1,
            const point& P
        );


public:

    //- Return types for classify
    enum proxType
    {
        NONE, 
        POINT,  // Close to point
        EDGE    // Close to edge
    };


    // Constructors

        //- Construct from three points
        inline triangle(const Point& a, const Point& b, const Point& c);

        //- Construct from Istream
        inline triangle(Istream&);


    // Member Functions

        // Access

            //- Return first vertex
            inline const Point& a() const;

            //- Return second vertex
            inline const Point& b() const;

            //- Return third vertex
            inline const Point& c() const;


        // Properties

            //- Return centre (centroid)
            inline Point centre() const;

            //- Return scalar magnitude
            inline scalar mag() const;

            //- Return vector normal
            inline vector normal() const;

            //- Return circum-centre
            inline vector circumCentre() const;

            //- Return circum-radius
            inline scalar circumRadius() const;

            //- Return quality: Ratio triangle and circum-circle area
            inline scalar quality() const;

            //- Return swept-volume
            inline scalar sweptVol(const triangle& t) const;

            //- Return point intersection with a ray.
            //  For a hit, the distance is signed. Positive number
            //  represents the point in front of triangle.
            //  In case of miss pointHit is set to nearest point
            //  on triangle and its distance to the distance between
            //  the original point and the plane intersection point
            inline pointHit ray
            (
                const point& p,
                const vector& q,
                const intersection::algorithm = intersection::FULL_RAY,
                const intersection::direction dir = intersection::VECTOR
            ) const;

            //- Fast intersection with a ray.
            //  For a hit, the pointHit.distance() is the line parameter t :
            //  intersection=p+t*q. Only defined for VISIBLE, FULL_RAY or
            //  HALF_RAY. tol increases the virtual size of the triangle
            // by a relative factor.
            inline pointHit intersection
            (
                const point& p,
                const vector& q,
                const intersection::algorithm alg,
                const scalar tol = 0.0
            ) const;

            //- Return nearest point to p on triangle
            inline pointHit nearestPoint
            (
                const point& p
            ) const;

            //- Classify point in triangle plane w.r.t. triangle edges.
            //  - inside (true returned)/outside (false returned)
            //  - near point (nearType=POINT, nearLabel=0, 1, 2)
            //  - near edge (nearType=EDGE, nearLabel=0, 1, 2)
            //    Note: edges are counted from starting
            //    vertex so e.g. edge 2 is from f[2] to f[0]
            //  tol is fraction to account for truncation error. Is only used
            //  when comparing normalized (0..1) numbers.
            bool classify
            (
                const point& p,
                const scalar tol,
                label& nearType,
                label& nearLabel
            ) const;


    // IOstream operators

        friend Istream& operator>> <Point, PointRef>
        (
            Istream&,
            triangle&
        );

        friend Ostream& operator<< <Point, PointRef>
        (
            Ostream&,
            const triangle&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triangleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
