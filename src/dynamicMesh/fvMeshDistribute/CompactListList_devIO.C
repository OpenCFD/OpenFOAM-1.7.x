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

#include "CompactListList_dev.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T, class Container>
Foam::CompactListList_dev<T, Container>::CompactListList_dev(Istream& is)
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, class Container>
Foam::Istream& Foam::operator>>(Istream& is, CompactListList_dev<T, Container>& lst)
{
    is  >> lst.offsets_ >> lst.m_;
    // Note: empty list gets output as two empty lists
    if (lst.offsets_.size() == 0)
    {
        lst.size_ = 0;
    }
    else
    {
        lst.size_ = lst.offsets_.size()-1;
    }
    return is;
}


template<class T, class Container>
Foam::Ostream& Foam::operator<<(Ostream& os, const CompactListList_dev<T, Container>& lst)
{
    os  << lst.offsets_ << lst.m_;
    return os;
}


// ************************************************************************* //
