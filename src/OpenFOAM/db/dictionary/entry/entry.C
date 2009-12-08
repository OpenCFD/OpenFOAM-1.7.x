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

#include "entry.H"
#include "dictionary.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::entry::entry(const keyType& keyword)
:
    keyword_(keyword)
{}


Foam::entry::entry(const entry& e)
:
    keyword_(e.keyword_)
{}


Foam::autoPtr<Foam::entry> Foam::entry::clone() const
{
    return clone(dictionary::null);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::entry::operator=(const entry& e)
{
    // check for assignment to self
    if (this == &e)
    {
        FatalErrorIn("entry::operator=(const entry&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    keyword_ = e.keyword_;
}


bool Foam::entry::operator==(const entry& e) const
{
    if (keyword_ != e.keyword_)
    {
        return false;
    }
    else
    {
        OStringStream oss1;
        oss1 << *this;

        OStringStream oss2;
        oss2 << e;

        return oss1.str() == oss2.str();
    }
}


bool Foam::entry::operator!=(const entry& e) const
{
    return !operator==(e);
}


// ************************************************************************* //
