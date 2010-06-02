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

\*---------------------------------------------------------------------------*/

#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if INT_MAX    != 2147483647
#    error "INT_MAX    != 2147483647"
#    error "The random number generator may not work!"
#endif

#ifdef USE_RANDOM
#   include <climits>
#else
#   include <cstdlib>
#endif


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// construct given seed
Random::Random(const label& seed)
{
    if (seed > 1)
    {
        Seed = seed;
    }
    else
    {
        Seed = 1;
    }

#   ifdef USE_RANDOM
        srandom((unsigned int)Seed);
#   else
        srand48(Seed);
#   endif

}


int Random::bit()
{
#   ifdef USE_RANDOM
    if (random() > INT_MAX/2)
#   else
    if (lrand48() > INT_MAX/2)
#   endif
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


scalar Random::scalar01()
{
#   ifdef USE_RANDOM
        return (scalar)random()/INT_MAX;
#   else
        return drand48();
#   endif
}


vector Random::vector01()
{
    vector rndVec;
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        rndVec.component(cmpt) = scalar01();
    }

    return rndVec;
}


sphericalTensor Random::sphericalTensor01()
{
    sphericalTensor rndTen;
    rndTen.ii() = scalar01();

    return rndTen;
}


symmTensor Random::symmTensor01()
{
    symmTensor rndTen;
    for (direction cmpt=0; cmpt<symmTensor::nComponents; cmpt++)
    {
        rndTen.component(cmpt) = scalar01();
    }

    return rndTen;
}


tensor Random::tensor01()
{
    tensor rndTen;
    for (direction cmpt=0; cmpt<tensor::nComponents; cmpt++)
    {
        rndTen.component(cmpt) = scalar01();
    }

    return rndTen;
}


label Random::integer(const label lower, const label upper)
{
#   ifdef USE_RANDOM
        return lower + (random() % (upper+1-lower));
#   else
        return lower + (lrand48() % (upper+1-lower));
#   endif
}


vector Random::position(const vector& start, const vector& end)
{
    vector rndVec(start);

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        rndVec.component(cmpt) +=
            scalar01()*(end.component(cmpt) - start.component(cmpt));
    }

    return rndVec;
}


void Random::randomise(scalar& s)
{
     s = scalar01();
}


void Random::randomise(vector& v)
{
    v = vector01();
}


void Random::randomise(sphericalTensor& st)
{
    st = sphericalTensor01();
}


void Random::randomise(symmTensor& st)
{
    st = symmTensor01();
}


void Random::randomise(tensor& t)
{
    t = tensor01();
}


// return a normal Gaussian randon number
// with zero mean and unity variance N(0, 1)

scalar Random::GaussNormal()
{
    static int iset = 0;
    static scalar gset;
    scalar fac, rsq, v1, v2;

    if (iset == 0)
    {
        do
        {
            v1 = 2.0*scalar01() - 1.0;
            v2 = 2.0*scalar01() - 1.0;
            rsq = v1*v1 + v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);

        fac = sqrt(-2.0 * log(rsq)/rsq);
        gset = v1*fac;
        iset = 1;

        return v2*fac;
    }
    else
    {
        iset = 0;

        return gset;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
