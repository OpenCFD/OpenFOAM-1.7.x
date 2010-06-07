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

Description

\*---------------------------------------------------------------------------*/

#include "triangle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline tetrahedron<Point, PointRef>::tetrahedron
(
    const Point& a,
    const Point& b,
    const Point& c,
    const Point& d
)
:
    a_(a),
    b_(b),
    c_(c),
    d_(d)
{}


template<class Point, class PointRef>
inline tetrahedron<Point, PointRef>::tetrahedron(Istream& is)
{
    // Read beginning of tetrahedron point pair
    is.readBegin("tetrahedron");

    is >> a_ >> b_ >> c_ >> d_;

    // Read end of tetrahedron point pair
    is.readEnd("tetrahedron");

    // Check state of Istream
    is.check("tetrahedron::tetrahedron(Istream& is)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Point, class PointRef>
inline const Point& tetrahedron<Point, PointRef>::a() const
{
    return a_;
}


template<class Point, class PointRef>
inline const Point& tetrahedron<Point, PointRef>::b() const
{
    return b_;
}


template<class Point, class PointRef>
inline const Point& tetrahedron<Point, PointRef>::c() const
{
    return c_;
}


template<class Point, class PointRef>
inline const Point& tetrahedron<Point, PointRef>::d() const
{
    return d_;
}


template<class Point, class PointRef>
inline vector tetrahedron<Point, PointRef>::Sa() const
{
    return triangle<Point, PointRef>(b_, c_, d_).normal();
}


template<class Point, class PointRef>
inline vector tetrahedron<Point, PointRef>::Sb() const
{
    return triangle<Point, PointRef>(a_, d_, c_).normal();
}


template<class Point, class PointRef>
inline vector tetrahedron<Point, PointRef>::Sc() const
{
    return triangle<Point, PointRef>(a_, b_, d_).normal();
}


template<class Point, class PointRef>
inline vector tetrahedron<Point, PointRef>::Sd() const
{
    return triangle<Point, PointRef>(a_, c_, b_).normal();
}


template<class Point, class PointRef>
inline scalar tetrahedron<Point, PointRef>::mag() const
{
    return (1.0/6.0)*(((b_ - a_) ^ (c_ - a_)) & (d_ - a_));
}


template<class Point, class PointRef>
inline vector tetrahedron<Point, PointRef>::circumCentre() const
{
    vector a = b_ - a_;
    vector b = c_ - a_;
    vector c = d_ - a_;

    scalar lamda = magSqr(c) - (a&c);
    scalar mu = magSqr(b) - (a&b);

    vector ba = b^a;
    vector ca = c^a;

    return a_ + 0.5*(a + (lamda*ba - mu*ca)/(c&ba));
}


template<class Point, class PointRef>
inline scalar tetrahedron<Point, PointRef>::circumRadius() const
{
    return Foam::mag(a_ - circumCentre());
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class point, class pointRef>
inline Istream& operator>>(Istream& is, tetrahedron<point, pointRef>& t)
{
    // Read beginning of tetrahedron point pair
    is.readBegin("tetrahedron");

    is >> t.a_ >> t.b_ >> t.c_ >> t.d_;

    // Read end of tetrahedron point pair
    is.readEnd("tetrahedron");

    // Check state of Ostream
    is.check("Istream& operator>>(Istream&, tetrahedron&)");

    return is;
}


template<class Point, class PointRef>
inline Ostream& operator<<(Ostream& os, const tetrahedron<Point, PointRef>& t)
{
    os  << nl
        << token::BEGIN_LIST
        << t.a_ << token::SPACE << t.b_
        << token::SPACE << t.c_ << token::SPACE << t.d_
        << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
