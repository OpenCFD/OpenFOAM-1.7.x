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

#ifndef HashSet_C
#define HashSet_C

#include "HashSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Key, class Hash>
template<class AnyType, class AnyHash>
Foam::HashSet<Key, Hash>::HashSet
(
    const HashTable<AnyType, Key, AnyHash>& h
)
:
    HashTable<nil, Key, Hash>(h.size())
{
    for
    (
        typename HashTable<AnyType, Key, AnyHash>::const_iterator
        cit = h.cbegin();
        cit != h.cend();
        ++cit
    )
    {
        insert(cit.key());
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Key, class Hash>
inline bool Foam::HashSet<Key, Hash>::operator[](const Key& key) const
{
    return found(key);
}


template<class Key, class Hash>
bool Foam::HashSet<Key, Hash>::operator==(const HashSet<Key, Hash>& rhs) const
{
    // Are all lhs elements in rhs?
    for (const_iterator iter = this->cbegin(); iter != this->cend(); ++iter)
    {
        if (!rhs.found(iter.key()))
        {
            return false;
        }
    }

    // Are all rhs elements in lhs?
    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        if (!found(iter.key()))
        {
            return false;
        }
    }

    return true;
}


template<class Key, class Hash>
bool Foam::HashSet<Key, Hash>::operator!=(const HashSet<Key, Hash>& rhs) const
{
    return !(operator==(rhs));
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::operator|=(const HashSet<Key, Hash>& rhs)
{
    // Add rhs elements into lhs
    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        insert(iter.key());
    }
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::operator&=(const HashSet<Key, Hash>& rhs)
{
    // Remove elements not also found in rhs
    for (iterator iter = this->begin(); iter != this->end(); ++iter)
    {
        if (!rhs.found(iter.key()))
        {
            erase(iter);
        }
    }
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::operator^=(const HashSet<Key, Hash>& rhs)
{
    // Add missed rhs elements, remove duplicate elements
    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        if (found(iter.key()))
        {
            erase(iter.key());
        }
        else
        {
            insert(iter.key());
        }
    }
}


// same as HashTable::erase()
template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::operator-=(const HashSet<Key, Hash>& rhs)
{
    // Remove rhs elements from lhs
    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        erase(iter.key());
    }
}


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

template<class Key, class Hash>
Foam::HashSet<Key, Hash>
Foam::operator|
(
    const HashSet<Key, Hash>& hash1,
    const HashSet<Key, Hash>& hash2
)
{
    HashSet<Key, Hash> out(hash1);
    out |= hash2;
    return out;
}


template<class Key, class Hash>
Foam::HashSet<Key, Hash>
Foam::operator&
(
    const HashSet<Key, Hash>& hash1,
    const HashSet<Key, Hash>& hash2
)
{
    HashSet<Key, Hash> out(hash1);
    out &= hash2;
    return out;
}


template<class Key, class Hash>
Foam::HashSet<Key, Hash>
Foam::operator^
(
    const HashSet<Key, Hash>& hash1,
    const HashSet<Key, Hash>& hash2
)
{
    HashSet<Key, Hash> out(hash1);
    out ^= hash2;
    return out;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
