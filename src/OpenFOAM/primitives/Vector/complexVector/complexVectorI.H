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

Description
    complexVector specific part of 3D complexVector obtained from
    generic Vector.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline complexVector operator*(const complex& v1, const complexVector& v2)
{
    return complexVector
    (
        v1*v2.x(),
        v1*v2.y(),
        v1*v2.z()
    );
}


inline complexVector operator*(const complexVector& v2, const complex& v1)
{
    return complexVector
    (
        v1*v2.x(),
        v1*v2.y(),
        v1*v2.z()
    );
}


inline complexVector operator/(const complexVector& v1, const complex& v2)
{
    return complexVector
    (
        v1.x()/v2,
        v1.y()/v2,
        v1.z()/v2
    );
}


inline complexVector operator/(const complex& v1, const complexVector& v2)
{
    return complexVector
    (
        v1/v2.x(),
        v1/v2.y(),
        v1/v2.z()
    );
}


// complexVector dot product

inline complex operator&(const complexVector& v1, const complexVector& v2)
{
    return complex
    (
        v1.x()*v2.x().conjugate()
      + v1.y()*v2.y().conjugate()
      + v1.z()*v2.z().conjugate()
    );
}


// complexVector cross product

inline complexVector operator^(const complexVector& v1, const complexVector& v2)
{
    return complexVector
    (
        (v1.y()*v2.z() - v1.z()*v2.y()),
        (v1.z()*v2.x() - v1.x()*v2.z()),
        (v1.x()*v2.y() - v1.y()*v2.x())
    );
}


// complexVector cross product

inline complexVector operator^(const vector& v1, const complexVector& v2)
{
    return complexVector
    (
        (v1.y()*v2.z() - v1.z()*v2.y()),
        (v1.z()*v2.x() - v1.x()*v2.z()),
        (v1.x()*v2.y() - v1.y()*v2.x())
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
