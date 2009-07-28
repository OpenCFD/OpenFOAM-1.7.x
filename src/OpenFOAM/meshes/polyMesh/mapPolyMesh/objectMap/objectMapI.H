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

Description

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline objectMap::objectMap()
:
    index_(-1),
    masterObjects_(0)
{}


inline objectMap::objectMap(const label index, const labelList& master)
:
    index_(index),
    masterObjects_(master)
{}


inline objectMap::objectMap(Istream& is)
{
    // Read beginning of objectMap
    is.readBegin("objectMap");

    is >> index_ >> static_cast<labelList&>(masterObjects_);

    // Read master of objectMap
    is.readEnd("objectMap");

    // Check state of Istream
    is.check("objectMap::objectMap(Istream&)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

label& objectMap::index()
{
    return index_;
}


inline label objectMap::index() const
{
    return index_;
}


inline labelList& objectMap::masterObjects()
{
    return masterObjects_;
}


inline const labelList& objectMap::masterObjects() const
{
    return masterObjects_;
}


// * * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * //

inline bool operator==(const objectMap& a, const objectMap& b)
{
    return
    (
        (a.index_ == b.index_) && (a.masterObjects_ == b.masterObjects_)
    );
}


inline bool operator!=(const objectMap& a, const objectMap& b)
{
    return (!(a == b));
}


// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

inline Ostream& operator<<(Ostream& os, const objectMap& a)
{
    os  << token::BEGIN_LIST
        << a.index_ << token::SPACE
        << a.masterObjects_
        << token::END_LIST;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const objectMap&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // Master namespace Foam

// ************************************************************************* //
