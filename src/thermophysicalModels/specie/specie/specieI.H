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

\*---------------------------------------------------------------------------*/

#include "specie.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Construct from components
inline specie::specie
(
    const word& name,
    const scalar nMoles,
    const scalar molWeight
)
:
    name_(name),
    nMoles_(nMoles),
    molWeight_(molWeight)
{}


// Construct from components without name
inline specie::specie
(
    const scalar nMoles,
    const scalar molWeight
)
:
    //name_(),
    nMoles_(nMoles),
    molWeight_(molWeight)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct as copy
inline specie::specie(const specie& st)
:
    name_(st.name_),
    nMoles_(st.nMoles_),
    molWeight_(st.molWeight_)
{}


// Construct as named copy
inline specie::specie(const word& name, const specie& st)
:
    name_(name),
    nMoles_(st.nMoles_),
    molWeight_(st.molWeight_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Molecular weight [kg/kmol]
inline scalar specie::W() const
{
    return molWeight_;
}

//- No of moles of this species in mixture
inline scalar specie::nMoles() const
{
    return nMoles_;
}

//- Gas constant [J/(kg K)]
inline scalar specie::R() const
{
    return RR/molWeight_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void specie::operator=(const specie& st)
{
  //name_ = st.name_;
    nMoles_ = st.nMoles_;
    molWeight_ = st.molWeight_;
}


inline void specie::operator+=(const specie& st)
{
    scalar sumNmoles_ = max(nMoles_ + st.nMoles_, SMALL);

    molWeight_ =
        nMoles_/sumNmoles_*molWeight_
      + st.nMoles_/sumNmoles_*st.molWeight_;

    nMoles_ = sumNmoles_;
}


inline void specie::operator-=(const specie& st)
{
    scalar diffnMoles_ = nMoles_ - st.nMoles_;
    if (mag(diffnMoles_) < SMALL)
    {
        diffnMoles_ = SMALL;
    }

    molWeight_ =
        nMoles_/diffnMoles_*molWeight_
      - st.nMoles_/diffnMoles_*st.molWeight_;

    nMoles_ = diffnMoles_;
}


inline void specie::operator*=(const scalar s)
{
    nMoles_ *= s;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline specie operator+(const specie& st1, const specie& st2)
{
    scalar sumNmoles_ = max(st1.nMoles_ + st2.nMoles_, SMALL);

    return specie
    (
        sumNmoles_,
        st1.nMoles_/sumNmoles_*st1.molWeight_
      + st2.nMoles_/sumNmoles_*st2.molWeight_
    );
}


inline specie operator-(const specie& st1, const specie& st2)
{
    scalar diffNmoles_ = st1.nMoles_ - st2.nMoles_;
    if (mag(diffNmoles_) < SMALL)
    {
        diffNmoles_ = SMALL;
    }

    return specie
    (
        diffNmoles_,
        st1.nMoles_/diffNmoles_*st1.molWeight_
      - st2.nMoles_/diffNmoles_*st2.molWeight_
    );
}


inline specie operator*(const scalar s, const specie& st)
{
    return specie
    (
        s*st.nMoles_,
        st.molWeight_
    );
}


inline specie operator==(const specie& st1, const specie& st2)
{
    return st2 - st1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
