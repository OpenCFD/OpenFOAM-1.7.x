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

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "pointHit.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Point, class PointRef>
pointHit triangle<Point, PointRef>::nearestPoint
(
    const Point& baseVertex,
    const vector& E0,
    const vector& E1,
    const point& P
)
{
    // Distance vector
    const vector D(baseVertex - P);

    // Some geometrical factors
    const scalar a = E0 & E0;
    const scalar b = E0 & E1;
    const scalar c = E1 & E1;

    // Precalculate distance factors
    const scalar d = E0 & D;
    const scalar e = E1 & D;
    const scalar f = D & D;

    // Do classification
    const scalar det = a*c - b*b;
    scalar s = b*e - c*d;
    scalar t = b*d - a*e;

    bool inside = false;

    if (s+t < det)
    {
        if (s < 0)
        {
            if (t < 0)
            {
                // Region 4
                if (e > 0)
                {
                    // min on edge t = 0
                    t = 0;
                    s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
                }
                else
                {
                    // min on edge s=0
                    s = 0;
                    t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
                }
            }
            else
            {
                // Region 3. Min on edge s = 0
                s = 0;
                t = (e >= 0 ? 0 : (-e >= c ? 1 : -e/c));
            }
        }
        else if (t < 0)
        {
            // Region 5
            t = 0;
            s = (d >= 0 ? 0 : (-d >= a ? 1 : -d/a));
        }
        else
        {
            // Region 0
            const scalar invDet = 1/det;
            s *= invDet;
            t *= invDet;

            inside = true;
        }
    }
    else
    {
        if (s < 0)
        {
            // Region 2
            const scalar tmp0 = b + d;
            const scalar tmp1 = c + e;
            if (tmp1 > tmp0)
            {
                // min on edge s+t=1
                const scalar numer = tmp1 - tmp0;
                const scalar denom = a-2*b+c;
                s = (numer >= denom ? 1 : numer/denom);
                t = 1 - s;
            }
            else
            {
                // min on edge s=0
                s = 0;
                t = (tmp1 <= 0 ? 1 : (e >= 0 ? 0 : - e/c));
            }
        }
        else if (t < 0)
        {
            // Region 6
            const scalar tmp0 = b + d;
            const scalar tmp1 = c + e;
            if (tmp1 > tmp0)
            {
                // min on edge s+t=1
                const scalar numer = tmp1 - tmp0;
                const scalar denom = a-2*b+c;
                s = (numer >= denom ? 1 : numer/denom);
                t = 1 - s;
            }
            else
            {
                // min on edge t=0
                t = 0;
                s = (tmp1 <= 0 ? 1 : (d >= 0 ? 0 : - d/a));
            }
        }
        else
        {
            // Region 1
            const scalar numer = c+e-(b+d);
            if (numer <= 0)
            {
                s = 0;
            }
            else
            {
                const scalar denom = a-2*b+c;
                s = (numer >= denom ? 1 : numer/denom);
            }
        }

        t = 1 - s;
    }

    // Calculate distance.
    // Note: Foam::mag used since truncation error causes negative distances
    // with points very close to one of the triangle vertices.
    // (Up to -2.77556e-14 seen). Could use +SMALL but that not large enough.

    return pointHit
    (
        inside,
        baseVertex + s*E0 + t*E1,
        Foam::sqrt
        (
            Foam::mag(a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f)
        ),
        !inside
    );
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline triangle<Point, PointRef>::triangle
(
    const Point& a,
    const Point& b,
    const Point& c
)
:
    a_(a),
    b_(b),
    c_(c)
{}


template<class Point, class PointRef>
inline triangle<Point, PointRef>::triangle(Istream& is)
{
    // Read beginning of triangle point pair
    is.readBegin("triangle");

    is >> a_ >> b_ >> c_;

    // Read end of triangle point pair
    is.readEnd("triangle");

    // Check state of Istream
    is.check("triangle::triangle(Istream& is)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline const Point& triangle<Point, PointRef>::a() const
{
    return a_;
}

template<class Point, class PointRef>
inline const Point& triangle<Point, PointRef>::b() const
{
    return b_;
}

template<class Point, class PointRef>
inline const Point& triangle<Point, PointRef>::c() const
{
    return c_;
}


template<class Point, class PointRef>
inline Point triangle<Point, PointRef>::centre() const
{
    return (1.0/3.0)*(a_ + b_ + c_);
}


template<class Point, class PointRef>
inline scalar triangle<Point, PointRef>::mag() const
{
    return ::Foam::mag(normal());
}


template<class Point, class PointRef>
inline vector triangle<Point, PointRef>::normal() const
{
    return 0.5*((b_ - a_)^(c_ - a_));
}


template<class Point, class PointRef>
inline vector triangle<Point, PointRef>::circumCentre() const
{
    scalar d1 = (c_ - a_)&(b_ - a_);
    scalar d2 = -(c_ - b_)&(b_ - a_);
    scalar d3 = (c_ - a_)&(c_ - b_);

    scalar c1 = d2*d3;
    scalar c2 = d3*d1;
    scalar c3 = d1*d2;

    scalar c = c1 + c2 + c3;

    return
    (
        ((c2 + c3)*a_ + (c3 + c1)*b_ + (c1 + c2)*c_)/(2*c)
    );
}


template<class Point, class PointRef>
inline scalar triangle<Point, PointRef>::circumRadius() const
{
    scalar d1 = (c_ - a_) & (b_ - a_);
    scalar d2 = - (c_ - b_) & (b_ - a_);
    scalar d3 = (c_ - a_) & (c_ - b_);

    scalar denom = d2*d3 + d3*d1 + d1*d2;

    if (Foam::mag(denom) < VSMALL)
    {
        return GREAT;
    }
    else
    {
        scalar a = (d1 + d2)*(d2 + d3)*(d3 + d1) / denom;

        return 0.5*Foam::sqrt(min(GREAT, max(0, a)));
    }
}


template<class Point, class PointRef>
inline scalar triangle<Point, PointRef>::quality() const
{
    return
        mag()
      / (
            mathematicalConstant::pi
          * Foam::sqr(circumRadius())
          * 0.413497
          + VSMALL
        );
}


template<class Point, class PointRef>
inline scalar triangle<Point, PointRef>::sweptVol(const triangle& t) const
{
    return (1.0/12.0)*
    (
        ((t.a_ - a_) & ((b_ - a_)^(c_ - a_)))
      + ((t.b_ - b_) & ((c_ - b_)^(t.a_ - b_)))
      + ((c_ - t.c_) & ((t.b_ - t.c_)^(t.a_ - t.c_)))

      + ((t.a_ - a_) & ((b_ - a_)^(c_ - a_)))
      + ((b_ - t.b_) & ((t.a_ - t.b_)^(t.c_ - t.b_)))
      + ((c_ - t.c_) & ((b_ - t.c_)^(t.a_ - t.c_)))
    );
}


template<class Point, class PointRef>
inline pointHit triangle<Point, PointRef>::ray
(
    const point& p,
    const vector& q,
    const intersection::algorithm alg,
    const intersection::direction dir
) const
{
    // Express triangle in terms of baseVertex (point a_) and
    // two edge vectors
    vector E0 = b_ - a_;
    vector E1 = c_ - a_;

    // Initialize intersection to miss.
    pointHit inter(p);

    vector n(0.5*(E0 ^ E1));
    scalar area = Foam::mag(n);

    if (area < VSMALL)
    {
        // Ineligible miss.
        inter.setMiss(false);

        // The miss point is the nearest point on the triangle. Take any one.
        inter.setPoint(a_);

        // The distance to the miss is the distance between the
        // original point and plane of intersection. No normal so use
        // distance from p to a_
        inter.setDistance(Foam::mag(a_ - p));

        return inter;
    }

    vector q1 = q/Foam::mag(q);

    if (dir == intersection::CONTACT_SPHERE)
    {
        n /= area;

        return ray(p, q1 - n, alg, intersection::VECTOR);
    }

    point pInter;
    bool hit;
    {
        // Reuse the fast ray intersection routine below in FULL_RAY
        // mode since the original intersection routine has rounding problems.
        pointHit fastInter = intersection(p, q1, intersection::FULL_RAY);
        pInter = fastInter.rawPoint();
        hit = fastInter.hit();
    }

    scalar dist = q1 & (pInter - p);

    const scalar planarPointTol =
        Foam::min
        (
            Foam::min
            (
                Foam::mag(E0),
                Foam::mag(E1)
            ),
            Foam::mag(c_ - b_)
        )*intersection::planarTol();

    bool eligible =
        alg == intersection::FULL_RAY
     || (alg == intersection::HALF_RAY && dist > -planarPointTol)
     || (
            alg == intersection::VISIBLE
         && ((q1 & normal()) < -VSMALL)
        );

    if (hit && eligible)
    {
        // Hit. Set distance to intersection.
        inter.setHit();
        inter.setPoint(pInter);
        inter.setDistance(dist);
    }
    else
    {
        // Miss or ineligible hit.
        inter.setMiss(eligible);

        // The miss point is the nearest point on the triangle
        inter.setPoint(nearestPoint(a_, E0, E1, p).rawPoint());

        // The distance to the miss is the distance between the
        // original point and plane of intersection
        inter.setDistance(Foam::mag(pInter - p));
    }


    return inter;
}


// From "Fast, Minimum Storage Ray/Triangle Intersection"
// Moeller/Trumbore.
template<class Point, class PointRef>
inline pointHit triangle<Point, PointRef>::intersection
(
    const point& orig,
    const vector& dir,
    const intersection::algorithm alg,
    const scalar tol
) const
{
    const vector edge1 = b_ - a_;
    const vector edge2 = c_ - a_;

    // begin calculating determinant - also used to calculate U parameter
    const vector pVec = dir ^ edge2;

    // if determinant is near zero, ray lies in plane of triangle
    const scalar det = edge1 & pVec;

    // Initialise to miss
    pointHit intersection(false, vector::zero, GREAT, false);

    if (alg == intersection::VISIBLE)
    {
        // Culling branch
        if (det < ROOTVSMALL)
        {
            // Ray on wrong side of triangle. Return miss
            return intersection;
        }
    }
    else if (alg == intersection::HALF_RAY || alg == intersection::FULL_RAY)
    {
        // Non-culling branch
        if (det > -ROOTVSMALL && det < ROOTVSMALL)
        {
            // Ray parallel to triangle. Return miss
            return intersection;
        }
    }

    const scalar inv_det = 1.0 / det;

    /* calculate distance from a_ to ray origin */
    const vector tVec = orig-a_;

    /* calculate U parameter and test bounds */
    const scalar u = (tVec & pVec)*inv_det;

    if (u < -tol || u > 1.0+tol)
    {
        // return miss
        return intersection;
    }

    /* prepare to test V parameter */
    const vector qVec = tVec ^ edge1;

    /* calculate V parameter and test bounds */
    const scalar v = (dir & qVec) * inv_det;

    if (v < -tol || u + v > 1.0+tol)
    {
        // return miss
        return intersection;
    }

    /* calculate t, scale parameters, ray intersects triangle */
    const scalar t = (edge2 & qVec) * inv_det;

    if (alg == intersection::HALF_RAY && t < -tol)
    {
        // Wrong side of orig. Return miss
        return intersection;
    }

    intersection.setHit();
    intersection.setPoint(a_ + u*edge1 + v*edge2);
    intersection.setDistance(t);

    return intersection;
}


template<class Point, class PointRef>
inline pointHit triangle<Point, PointRef>::nearestPoint
(
    const point& p
) const
{
    // Express triangle in terms of baseVertex (point a_) and
    // two edge vectors
    vector E0 = b_ - a_;
    vector E1 = c_ - a_;

    return nearestPoint(a_, E0, E1, p);
}


template<class Point, class PointRef>
inline bool triangle<Point, PointRef>::classify
(
    const point& p,
    const scalar tol,
    label& nearType,
    label& nearLabel
) const
{
    const vector E0 = b_ - a_;
    const vector E1 = c_ - a_;
    const vector n = 0.5*(E0 ^ E1);

    // Get largest component of normal
    scalar magX = Foam::mag(n.x());
    scalar magY = Foam::mag(n.y());
    scalar magZ = Foam::mag(n.z());

    label i0 = -1;
    if ((magX >= magY) && (magX >= magZ))
    {
        i0 = 0;
    }
    else if ((magY >= magX) && (magY >= magZ))
    {
        i0 = 1;
    }
    else
    {
        i0 = 2;
    }

    // Get other components
    label i1 = (i0 + 1) % 3;
    label i2 = (i1 + 1) % 3;


    scalar u1 = E0[i1];
    scalar v1 = E0[i2];

    scalar u2 = E1[i1];
    scalar v2 = E1[i2];

    scalar det = v2*u1 - u2*v1;

    scalar u0 = p[i1] - a_[i1];
    scalar v0 = p[i2] - a_[i2];

    scalar alpha = 0;
    scalar beta = 0;

    bool hit = false;

    if (Foam::mag(u1) < ROOTVSMALL)
    {
        beta = u0/u2;

        alpha = (v0 - beta*v2)/v1;

        hit = ((alpha >= 0) && ((alpha + beta) <= 1));
    }
    else
    {
        beta = (v0*u1 - u0*v1)/det;

        alpha = (u0 - beta*u2)/u1;

        hit = ((alpha >= 0) && ((alpha + beta) <= 1));
    }

    //
    // Now alpha, beta are the coordinates in the triangle local coordinate
    // system E0, E1
    //

    //Pout<< "alpha:" << alpha << endl;
    //Pout<< "beta:" << beta << endl;
    //Pout<< "hit:" << hit << endl;
    //Pout<< "tol:" << tol << endl;

    if (hit)
    {
        // alpha,beta might get negative due to precision errors
        alpha = max(0.0, min(1.0, alpha));
        beta = max(0.0, min(1.0, beta));
    }

    nearType = NONE;
    nearLabel = -1;

    if (Foam::mag(alpha+beta-1) <= tol)
    {
        // On edge between vert 1-2 (=edge 1)

        if (Foam::mag(alpha) <= tol)
        {
            nearType = POINT;
            nearLabel = 2;
        }
        else if (Foam::mag(beta) <= tol)
        {
            nearType = POINT;
            nearLabel = 1;
        }
        else if ((alpha >= 0) && (alpha <= 1) && (beta >= 0) && (beta <= 1))
        {
            nearType = EDGE;
            nearLabel = 1;
        }
    }
    else if (Foam::mag(alpha) <= tol)
    {
        // On edge between vert 2-0 (=edge 2)

        if (Foam::mag(beta) <= tol)
        {
            nearType = POINT;
            nearLabel = 0;
        }
        else if (Foam::mag(beta-1) <= tol)
        {
            nearType = POINT;
            nearLabel = 2;
        }
        else if ((beta >= 0) && (beta <= 1))
        {
            nearType = EDGE;
            nearLabel = 2;
        }
    }
    else if (Foam::mag(beta) <= tol)
    {
        // On edge between vert 0-1 (= edge 0)

        if (Foam::mag(alpha) <= tol)
        {
            nearType = POINT;
            nearLabel = 0;
        }
        else if (Foam::mag(alpha-1) <= tol)
        {
            nearType = POINT;
            nearLabel = 1;
        }
        else if ((alpha >= 0) && (alpha <= 1))
        {
            nearType = EDGE;
            nearLabel = 0;
        }
    }

    return hit;
}





// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class point, class pointRef>
inline Istream& operator>>(Istream& is, triangle<point, pointRef>& t)
{
    // Read beginning of triangle point pair
    is.readBegin("triangle");

    is >> t.a_ >> t.b_ >> t.c_;

    // Read end of triangle point pair
    is.readEnd("triangle");

    // Check state of Ostream
    is.check("Istream& operator>>(Istream&, triangle&)");

    return is;
}


template<class Point, class PointRef>
inline Ostream& operator<<(Ostream& os, const triangle<Point, PointRef>& t)
{
    os  << nl
        << token::BEGIN_LIST
        << t.a_ << token::SPACE << t.b_ << token::SPACE << t.c_
        << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

