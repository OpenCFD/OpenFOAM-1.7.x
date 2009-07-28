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

Typedef
    Foam::floatScalar

Description
    single floating point number identical to float

SourceFiles
    floatScalar.C

\*---------------------------------------------------------------------------*/

#ifndef floatScalar_H
#define floatScalar_H

#include "label.H"
#include "word.H"

#include <cmath>

#ifdef ibm
    float lgamma(float);
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define floatScalar as a float

namespace Foam
{
    typedef float floatScalar;
}

/*
// Define floatScalar as a float

namespace Foam
{
    typedef float floatScalar;

    // Largest and smallest floatScalar values allowed in certain parts of the code
    // (6 is the number of significant figures in an
    //  IEEE single precision number.  See limits.h or float.h)
    static const floatScalar GREAT = 1.0e+6;
    static const floatScalar VGREAT = 1.0e+37;
    static const floatScalar SMALL = 1.0e-6;
    static const floatScalar VSMALL = 1.0e-37;
}
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pTraits.H"
#include "products.H"
#include "direction.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// template specialisation for pTraits<floatScalar>
template<>
class pTraits<floatScalar>
{
    floatScalar p_;

public:

    //- Component type
    typedef floatScalar cmptType;

    // Member constants

        enum
        {
            dim = 3,         // Dimensionality of space
            rank = 0,        // Rank od floatScalar is 0
            nComponents = 1  // Number of components in floatScalar is 1
        };

    // Static data members

        static const char* const typeName;
        static const char* componentNames[];
        static const floatScalar zero;
        static const floatScalar one;

    // Constructors

        //- Construct from Istream
        pTraits(Istream& is);

    // Member Functions

        operator floatScalar() const
        {
            return p_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//template<class Cmpt>
//class typeOfRank<Cmpt, 0>
//{
//public:
//
//    typedef Cmpt type;
//};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Return a string representation of a floatScalar
word name(const floatScalar s)
{
    return name(scalar(s));
}

#define MAXMINPOW(retType, type1, type2)          \
                                                  \
MAXMIN(retType, type1, type2)                     \
                                                  \
inline float pow(const type1 s, const type2 e)   \
{                                                 \
    return ::pow(float(s), float(e));           \
}


//MAXMINPOW(float, float, float)
//MAXMINPOW(float, float, int)
//MAXMINPOW(float, int, float)
//MAXMINPOW(float, float, long)
//MAXMINPOW(float, long, float)
//MAXMINPOW(float, float, int)
//MAXMINPOW(float, int, float)
//MAXMINPOW(float, float, long)
//MAXMINPOW(float, long, float)

#undef MAXMINPOW


inline floatScalar mag(const floatScalar s)
{
    return ::fabs(s);
}

inline floatScalar sign(const floatScalar s)
{
    return (s >= 0)? 1: -1;
}

inline floatScalar pos(const floatScalar s)
{
    return (s >= 0)? 1: 0;
}

inline floatScalar neg(const floatScalar s)
{
    return (s < 0)? 1: 0;
}

inline floatScalar limit(const floatScalar s1, const floatScalar s2)
{
    return (mag(s1) < mag(s2))? s1: 0.0;
}

inline floatScalar magSqr(const floatScalar s)
{
    return s*s;
}

inline floatScalar sqr(const floatScalar s)
{
    return s*s;
}

inline floatScalar pow3(const floatScalar s)
{
    return s*s*s;
}

inline floatScalar pow4(const floatScalar s)
{
    return sqr(sqr(s));
}

inline floatScalar cmptAv(const floatScalar s)
{
    return s;
}

inline floatScalar cmptMag(const floatScalar s)
{
    return mag(s);
}

inline floatScalar scale(const floatScalar s, const floatScalar d)
{
    return s*d;
}


#define transFunc(func)            \
inline floatScalar func(const floatScalar s) \
{                                  \
    return ::func(s);              \
}

// Standard C++ transcendental functions
transFunc(sqrt)
transFunc(exp)
transFunc(log)
transFunc(log10)
transFunc(sin)
transFunc(cos)
transFunc(tan)
transFunc(asin)
transFunc(acos)
transFunc(atan)
transFunc(sinh)
transFunc(cosh)
transFunc(tanh)
transFunc(asinh)
transFunc(acosh)
transFunc(atanh)

// Standard ANSI-C (but not in <cmath>) transcendental functions
transFunc(erf)
transFunc(erfc)
transFunc(lgamma)
transFunc(j0)
transFunc(j1)
transFunc(y0)
transFunc(y1)

#undef transFunc

// Stabilisation around zero for division
inline floatScalar stabilise(const floatScalar s, const floatScalar small)
{
    if (s >= 0)
    {
        return s + small;
    }
    else
    {
        return s - small;
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

//floatScalar readScalar(Istream& is);
Istream& operator>>(Istream&, floatScalar&);
Ostream& operator<<(Ostream&, const floatScalar);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
