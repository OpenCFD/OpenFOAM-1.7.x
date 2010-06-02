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

InClass
    Foam::Pstream

Description
    Combination-Reduction operation for a parallel run.

    The information from all nodes is collected on the master node,
    combined using the given combination function and the result is
    broadcast to all nodes

\*---------------------------------------------------------------------------*/

#ifndef ops_H
#define ops_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define EqOp(opName, op)                                                    \
                                                                            \
template<class T1, class T2>                                                \
class opName##Op2                                                           \
{                                                                           \
public:                                                                     \
                                                                            \
    void operator()(T1& x, const T2& y) const                               \
    {                                                                       \
        op;                                                                 \
    }                                                                       \
};                                                                          \
                                                                            \
template<class T>                                                           \
class opName##Op                                                            \
{                                                                           \
public:                                                                     \
                                                                            \
    void operator()(T& x, const T& y) const                                 \
    {                                                                       \
        op;                                                                 \
    }                                                                       \
};

EqOp(eq, x = y)
EqOp(plusEq, x += y)
EqOp(minusEq, x -= y)
EqOp(multiplyEq, x *= y)
EqOp(divideEq, x /= y)
EqOp(eqMag, x = mag(y))
EqOp(plusEqMagSqr, x += magSqr(y))
EqOp(maxEq, x = max(x, y))
EqOp(minEq, x = min(x, y))
EqOp(andEq, x = (x && y))
EqOp(orEq, x  = (x || y))

EqOp(eqMinus, x = -y)

#undef EqOp


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define Op(opName, op)                                                      \
                                                                            \
template<class T, class T1, class T2>                                       \
class opName##Op3                                                           \
{                                                                           \
public:                                                                     \
                                                                            \
    T operator()(const T1& x, const T2& y) const                            \
    {                                                                       \
        return op;                                                          \
    }                                                                       \
};                                                                          \
                                                                            \
template<class T1, class T2>                                                \
class opName##Op2                                                           \
{                                                                           \
public:                                                                     \
                                                                            \
    T1 operator()(const T1& x, const T2& y) const                           \
    {                                                                       \
        return op;                                                          \
    }                                                                       \
};                                                                          \
                                                                            \
template<class T>                                                           \
class opName##Op                                                            \
{                                                                           \
public:                                                                     \
                                                                            \
    T operator()(const T& x, const T& y) const                              \
    {                                                                       \
        return op;                                                          \
    }                                                                       \
};

Op(sum, x + y)

Op(plus, x + y)
Op(minus, x - y)
Op(multiply, x * y)
Op(divide, x / y)
Op(cmptMultiply, cmptMultiply(x, y))
Op(cmptDivide, cmptDivide(x, y))
Op(stabilise, stabilise(x, y))
Op(max, max(x, y))
Op(min, min(x, y))
Op(minMod, minMod(x, y))
Op(and, x && y)
Op(or, x || y)
Op(eqEq, x == y)

#undef Op


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
