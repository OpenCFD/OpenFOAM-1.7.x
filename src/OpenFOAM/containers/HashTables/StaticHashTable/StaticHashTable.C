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

#ifndef StaticHashTable_C
#define StaticHashTable_C

#include "StaticHashTable.H"
#include "List.H"
#include "IOstreams.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::label Foam::StaticHashTable<T, Key, Hash>::canonicalSize(const label size)
{
    if (size < 1)
    {
        return 0;
    }

    // enforce power of two
    unsigned int goodSize = size;

    if (goodSize & (goodSize - 1))
    {
        // brute-force is fast enough
        goodSize = 1;
        while (goodSize < unsigned(size))
        {
            goodSize <<= 1;
        }
    }

    return goodSize;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given initial table size
template<class T, class Key, class Hash>
Foam::StaticHashTable<T, Key, Hash>::StaticHashTable(const label size)
:
    StaticHashTableName(),
    keys_(canonicalSize(size)),
    objects_(keys_.size()),
    nElmts_(0),
    endIter_(*this, keys_.size(), 0),
    endConstIter_(*this, keys_.size(), 0)
{
    if (size < 1)
    {
        FatalErrorIn
        (
            "StaticHashTable<T, Key, Hash>::StaticHashTable(const label size)"
        )   << "Illegal size " << size << " for StaticHashTable."
            << " Minimum size is 1" << abort(FatalError);
    }
}


// Construct as copy
template<class T, class Key, class Hash>
Foam::StaticHashTable<T, Key, Hash>::StaticHashTable
(
    const StaticHashTable<T, Key, Hash>& ht
)
:
    StaticHashTableName(),
    keys_(ht.keys_),
    objects_(ht.objects_),
    nElmts_(ht.nElmts_),
    endIter_(*this, keys_.size(), 0),
    endConstIter_(*this, keys_.size(), 0)
{}



template<class T, class Key, class Hash>
Foam::StaticHashTable<T, Key, Hash>::StaticHashTable
(
    const Xfer< StaticHashTable<T, Key, Hash> >& ht
)
:
    StaticHashTableName(),
    keys_(0),
    objects_(0),
    nElmts_(0),
    endIter_(*this, 0, 0),
    endConstIter_(*this, 0, 0)
{
    transfer(ht());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
Foam::StaticHashTable<T, Key, Hash>::~StaticHashTable()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
bool Foam::StaticHashTable<T, Key, Hash>::found(const Key& key) const
{
    if (nElmts_)
    {
        const label hashIdx = hashKeyIndex(key);
        const List<Key>& localKeys = keys_[hashIdx];

        forAll(localKeys, elemIdx)
        {
            if (key == localKeys[elemIdx])
            {
                return true;
            }
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info<< "StaticHashTable<T, Key, Hash>::found(const Key&) : "
            << "Entry " << key << " not found in hash table\n";
    }
#   endif

    return false;
}


template<class T, class Key, class Hash>
typename Foam::StaticHashTable<T, Key, Hash>::iterator
Foam::StaticHashTable<T, Key, Hash>::find
(
    const Key& key
)
{
    if (nElmts_)
    {
        const label hashIdx = hashKeyIndex(key);
        const List<Key>& localKeys = keys_[hashIdx];

        forAll(localKeys, elemIdx)
        {
            if (key == localKeys[elemIdx])
            {
                return iterator(*this, hashIdx, elemIdx);
            }
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info<< "StaticHashTable<T, Key, Hash>::find(const Key&) : "
            << "Entry " << key << " not found in hash table\n";
    }
#   endif

    return end();
}


template<class T, class Key, class Hash>
typename Foam::StaticHashTable<T, Key, Hash>::const_iterator
Foam::StaticHashTable<T, Key, Hash>::find
(
    const Key& key
) const
{
    if (nElmts_)
    {
        const label hashIdx = hashKeyIndex(key);
        const List<Key>& localKeys = keys_[hashIdx];

        forAll(localKeys, elemIdx)
        {
            if (key == localKeys[elemIdx])
            {
                return const_iterator(*this, hashIdx, elemIdx);
            }
        }
    }

#   ifdef FULLDEBUG
    if (debug)
    {
        Info<< "StaticHashTable<T, Key, Hash>::find(const Key&) const : "
            << "Entry " << key << " not found in hash table\n";
    }
#   endif

    return cend();
}


// Return the table of contents
template<class T, class Key, class Hash>
Foam::List<Key> Foam::StaticHashTable<T, Key, Hash>::toc() const
{
    List<Key> tofc(nElmts_);
    label i = 0;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        tofc[i++] = iter.key();
    }

    return tofc;
}


template<class T, class Key, class Hash>
bool Foam::StaticHashTable<T, Key, Hash>::set
(
    const Key& key,
    const T& newEntry,
    const bool protect
)
{
    const label hashIdx = hashKeyIndex(key);
    List<Key>& localKeys = keys_[hashIdx];

    label existing = localKeys.size();
    forAll(localKeys, elemIdx)
    {
        if (key == localKeys[elemIdx])
        {
            existing = elemIdx;
            break;
        }
    }

    if (existing == localKeys.size())
    {
        // not found, append
        List<T>& localObjects = objects_[hashIdx];

        localKeys.setSize(existing+1);
        localObjects.setSize(existing+1);

        localKeys[existing] = key;
        localObjects[existing] = newEntry;

        nElmts_++;
    }
    else if (protect)
    {
        // found - but protected from overwriting
        // this corresponds to the STL 'insert' convention
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "StaticHashTable<T, Key, Hash>::set"
                "(const Key& key, T newEntry, true) : "
                "Cannot insert " << key << " already in hash table\n";
        }
#       endif
        return false;
    }
    else
    {
        // found - overwrite existing entry
        // this corresponds to the Perl convention
        objects_[hashIdx][existing] = newEntry;
    }

    return true;
}


template<class T, class Key, class Hash>
bool Foam::StaticHashTable<T, Key, Hash>::erase(const iterator& cit)
{
    if (cit != end())
    {
        List<Key>& localKeys = keys_[cit.hashIndex_];
        List<T>& localObjects = objects_[cit.hashIndex_];

        // Copy down
        for (label i = cit.elemIndex_+1; i < localKeys.size(); i++)
        {
            localKeys[i-1] = localKeys[i];
            localObjects[i-1] = localObjects[i];
        }
        localKeys.setSize(localKeys.size()-1);
        localObjects.setSize(localObjects.size()-1);

        // adjust iterator after erase
        iterator& it = const_cast<iterator&>(cit);

        it.elemIndex_--;
        if (it.elemIndex_ < 0)
        {
            // No previous element in the local list

            // Search back for previous non-zero table entry
            while (--it.hashIndex_ >= 0 && !objects_[it.hashIndex_].size())
            {}

            if (it.hashIndex_ >= 0)
            {
                // The last element in the local list
                it.elemIndex_ = objects_[it.hashIndex_].size() - 1;
            }
            else
            {
                // No previous found. Mark with special value which is
                // - not end()
                // - handled by operator++
                it.hashIndex_ = -1;
                it.elemIndex_ = 0;
            }
        }

        nElmts_--;

#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "StaticHashTable<T, Key, Hash>::erase(iterator&) : "
                << "hashedEntry removed.\n";
        }
#       endif

        return true;
    }
    else
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "StaticHashTable<T, Key, Hash>::erase(iterator&) : "
                << "cannot remove hashedEntry from hash table\n";
        }
#       endif

        return false;
    }
}


template<class T, class Key, class Hash>
bool Foam::StaticHashTable<T, Key, Hash>::erase(const Key& key)
{
    iterator it = find(key);

    if (it != end())
    {
        return erase(it);
    }
    else
    {
        return false;
    }
}


template<class T, class Key, class Hash>
Foam::label Foam::StaticHashTable<T, Key, Hash>::erase
(
    const StaticHashTable<T, Key, Hash>& rhs
)
{
    label count = 0;

    // Remove rhs elements from this table
    // NOTE: could optimize depending on which hash is smaller
    for (iterator iter = this->begin(); iter != this->end(); ++iter)
    {
        if (rhs.found(iter.key()) && erase(iter))
        {
            count++;
        }
    }

   return count;
}


template<class T, class Key, class Hash>
void Foam::StaticHashTable<T, Key, Hash>::resize(const label sz)
{
    label newSize = canonicalSize(sz);
    
    if (newSize == keys_.size())
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "StaticHashTable<T, Key, Hash>::resize(const label) : "
                << "new table size == old table size\n";
        }
#       endif

        return;
    }

    if (newSize < 1)
    {
        FatalErrorIn
        (
            "StaticHashTable<T, Key, Hash>::resize(const label)"
        )   << "Illegal size " << newSize << " for StaticHashTable."
            << " Minimum size is 1" << abort(FatalError);
    }


    StaticHashTable<T, Key, Hash> newTable(newSize);

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        newTable.insert(iter.key(), *iter);
    }

    transfer(newTable);

    // Adapt end() iterators
    endIter_.hashIndex_ = keys_.size();
    endConstIter_.hashIndex_ = keys_.size();
}


template<class T, class Key, class Hash>
void Foam::StaticHashTable<T, Key, Hash>::clear()
{
    forAll(keys_, hashIdx)
    {
        keys_[hashIdx].clear();
        objects_[hashIdx].clear();
    }

    nElmts_ = 0;
}


template<class T, class Key, class Hash>
void Foam::StaticHashTable<T, Key, Hash>::clearStorage()
{
    clear();
    resize(1);
}



template<class T, class Key, class Hash>
void Foam::StaticHashTable<T, Key, Hash>::transfer
(
    StaticHashTable<T, Key, Hash>& ht
)
{
    // Remove existing elements
    clear();

    // Copy data from ht
    keys_.transfer(ht.keys_);
    objects_.transfer(ht.objects_);

    nElmts_ = ht.nElmts_;
    ht.nElmts_ = 0;

    // Adapt end() iterators
    endIter_.hashIndex_ = keys_.size();
    endConstIter_.hashIndex_ = keys_.size();

    ht.endIter_.hashIndex_ = 0;
    ht.endConstIter_.hashIndex_ = 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
void Foam::StaticHashTable<T, Key, Hash>::operator=
(
    const StaticHashTable<T, Key, Hash>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "StaticHashTable<T, Key, Hash>::operator="
            "(const StaticHashTable<T, Key, Hash>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }


    // keys could be empty from a previous transfer()
    if (keys_.empty())
    {
        keys_.setSize(rhs.keys_.size());
        objects_.setSize(keys_.size());

        // Adapt end() iterators
        endIter_.hashIndex_ = keys_.size();
        endConstIter_.hashIndex_ = keys_.size();
    }
    else
    {
        clear();
        // keys_.size() does not change so neither does end() iterator.
    }


    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        insert(iter.key(), *iter);
    }
}

template<class T, class Key, class Hash>
bool Foam::StaticHashTable<T, Key, Hash>::operator==
(
    const StaticHashTable<T, Key, Hash>& rhs
) const
{
    // Are all my elements in rhs?
    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        const_iterator fnd = rhs.find(iter.key());

        if (fnd == rhs.cend() || fnd() != iter())
        {
            return false;
        }
    }

    // Are all rhs elements in me?
    for (const_iterator iter = rhs.cbegin(); iter != rhs.cend(); ++iter)
    {
        const_iterator fnd = find(iter.key());

        if (fnd == cend() || fnd() != iter())
        {
            return false;
        }
    }
    return true;
}


template<class T, class Key, class Hash>
bool Foam::StaticHashTable<T, Key, Hash>::operator!=
(
    const StaticHashTable<T, Key, Hash>& rhs
) const
{
    return !(operator==(rhs));
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#include "StaticHashTableIO.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
