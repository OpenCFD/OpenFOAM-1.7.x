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

#include "error.H"

// * * * * * * * * * * * * * Private Member Classes * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::hashedEntry::hashedEntry
(
    const Key& key,
    hashedEntry* next,
    const T& newEntry
)
:
    key_(key),
    next_(next),
    obj_(newEntry)
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::label
Foam::HashTable<T, Key, Hash>::hashKeyIndex(const Key& key) const
{
    // size is power of two - this is the modulus
    return Hash()(key) & (tableSize_ - 1);
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::label Foam::HashTable<T, Key, Hash>::size() const
{
    return nElmts_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::empty() const
{
    return !nElmts_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::insert
(
    const Key& key,
    const T& newEntry
)
{
    return set(key, newEntry, true);
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::set
(
    const Key& key,
    const T& newEntry
)
{
    return set(key, newEntry, false);
}


template<class T, class Key, class Hash>
inline Foam::Xfer<Foam::HashTable<T, Key, Hash> >
Foam::HashTable<T, Key, Hash>::xfer()
{
    return xferMove(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::operator[](const Key& key)
{
    iterator iter = find(key);

    if (iter == end())
    {
        FatalErrorIn("HashTable<T, Key, Hash>::operator[](const Key&)")
            << key << " not found in table.  Valid entries: "
            << toc()
            << exit(FatalError);
    }

    return *iter;
}


template<class T, class Key, class Hash>
inline const T& Foam::HashTable<T, Key, Hash>::operator[](const Key& key) const
{
    const_iterator iter = find(key);

    if (iter == cend())
    {
        FatalErrorIn("HashTable<T, Key, Hash>::operator[](const Key&) const")
            << key << " not found in table.  Valid entries: "
            << toc()
            << exit(FatalError);
    }

    return *iter;
}


template<class T, class Key, class Hash>
inline T& Foam::HashTable<T, Key, Hash>::operator()(const Key& key)
{
    iterator iter = find(key);

    if (iter == end())
    {
        insert(key, T());
        return *find(key);
    }
    else
    {
        return *iter;
    }
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::iterator::iterator
(
    HashTable<T, Key, Hash>& hashTbl,
    hashedEntry* elmt,
    label hashIndex
)
:
    hashTable_(hashTbl),
    elmtPtr_(elmt),
    hashIndex_(hashIndex)
{}


template<class T, class Key, class Hash>
inline void Foam::HashTable<T, Key, Hash>::iterator::operator=
(
    const iterator& iter
)
{
    elmtPtr_ = iter.elmtPtr_;
    hashIndex_ = iter.hashIndex_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::iterator::operator==
(
    const iterator& iter
) const
{
    return elmtPtr_ == iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::iterator::operator!=
(
    const iterator& iter
) const
{
    return elmtPtr_ != iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::iterator::operator==
(
    const const_iterator& iter
) const
{
    return elmtPtr_ == iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::iterator::operator!=
(
    const const_iterator& iter
) const
{
    return elmtPtr_ != iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline T&
Foam::HashTable<T, Key, Hash>::iterator::operator*()
{
    return elmtPtr_->obj_;
}


template<class T, class Key, class Hash>
inline T&
Foam::HashTable<T, Key, Hash>::iterator::operator()()
{
    return elmtPtr_->obj_;
}


template<class T, class Key, class Hash>
inline const T&
Foam::HashTable<T, Key, Hash>::iterator::operator*() const
{
    return elmtPtr_->obj_;
}


template<class T, class Key, class Hash>
inline const T&
Foam::HashTable<T, Key, Hash>::iterator::operator()() const
{
    return elmtPtr_->obj_;
}


template<class T, class Key, class Hash>
inline
typename Foam::HashTable<T, Key, Hash>::iterator&
Foam::HashTable<T, Key, Hash>::iterator::operator++()
{
    // Check for special value from erase. (sets hashIndex to -1)
    if (hashIndex_ >= 0)
    {
        // Do we have additional elements on the SLList?
        if (elmtPtr_ && elmtPtr_->next_)
        {
            elmtPtr_ = elmtPtr_->next_;
            return *this;
        }
    }

    // Step to the next table entry
    while
    (
        ++hashIndex_ < hashTable_.tableSize_
     && !(elmtPtr_ = hashTable_.table_[hashIndex_])
    )
    {}

    if (hashIndex_ == hashTable_.tableSize_)
    {
        // make end iterator
        elmtPtr_ = 0;
        hashIndex_ = 0;
    }
    return *this;
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::iterator
Foam::HashTable<T, Key, Hash>::iterator::operator++
(
    int
)
{
    iterator tmp = *this;
    ++*this;
    return tmp;
}


template<class T, class Key, class Hash>
inline
const Key& Foam::HashTable<T, Key, Hash>::iterator::key() const
{
    return elmtPtr_->key_;
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::iterator
Foam::HashTable<T, Key, Hash>::begin()
{
    label i = 0;

    if (nElmts_)
    {
        while (table_ && !table_[i] && ++i < tableSize_)
        {}
    }
    else
    {
        i = tableSize_;
    }

    if (i == tableSize_)
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "HashTable is empty\n";
        }
#       endif

        return HashTable<T, Key, Hash>::endIter_;
    }
    else
    {
        return iterator(*this, table_[i], i);
    }
}


template<class T, class Key, class Hash>
inline const typename Foam::HashTable<T, Key, Hash>::iterator&
Foam::HashTable<T, Key, Hash>::end()
{
    return HashTable<T, Key, Hash>::endIter_;
}


// * * * * * * * * * * * * * * * STL const_iterator * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::const_iterator::const_iterator
(
    const HashTable<T, Key, Hash>& hashTbl,
    const hashedEntry* elmt,
    label hashIndex
)
:
    hashTable_(hashTbl),
    elmtPtr_(elmt),
    hashIndex_(hashIndex)
{}


template<class T, class Key, class Hash>
inline Foam::HashTable<T, Key, Hash>::const_iterator::const_iterator
(
    const iterator& iter
)
:
    hashTable_(iter.hashTable_),
    elmtPtr_(iter.elmtPtr_),
    hashIndex_(iter.hashIndex_)
{}


template<class T, class Key, class Hash>
inline void Foam::HashTable<T, Key, Hash>::const_iterator::operator=
(
    const const_iterator& iter
)
{
    elmtPtr_ = iter.elmtPtr_;
    hashIndex_ = iter.hashIndex_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return elmtPtr_ == iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return elmtPtr_ != iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::const_iterator::operator==
(
    const iterator& iter
) const
{
    return elmtPtr_ == iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline bool Foam::HashTable<T, Key, Hash>::const_iterator::operator!=
(
    const iterator& iter
) const
{
    return elmtPtr_ != iter.elmtPtr_;
}


template<class T, class Key, class Hash>
inline const T&
Foam::HashTable<T, Key, Hash>::const_iterator::operator*() const
{
    return elmtPtr_->obj_;
}

template<class T, class Key, class Hash>
inline const T&
Foam::HashTable<T, Key, Hash>::const_iterator::operator()() const
{
    return elmtPtr_->obj_;
}


template<class T, class Key, class Hash>
inline
typename Foam::HashTable<T, Key, Hash>::const_iterator&
Foam::HashTable<T, Key, Hash>::const_iterator::operator++()
{
    if
    (
        !(elmtPtr_ = elmtPtr_->next_)
     && ++hashIndex_ < hashTable_.tableSize_
     && !(elmtPtr_ = hashTable_.table_[hashIndex_])
    )
    {
        while
        (
            ++hashIndex_ < hashTable_.tableSize_
         && !(elmtPtr_ = hashTable_.table_[hashIndex_])
        )
        {}
    }

    return *this;
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::const_iterator::operator++
(
    int
)
{
    const_iterator tmp = *this;
    ++*this;
    return tmp;
}


template<class T, class Key, class Hash>
inline
const Key& Foam::HashTable<T, Key, Hash>::const_iterator::key() const
{
    return elmtPtr_->key_;
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::cbegin() const
{
    label i = 0;

    if (nElmts_)
    {
        while (table_ && !table_[i] && ++i < tableSize_)
        {}
    }
    else
    {
        i = tableSize_;
    }

    if (i == tableSize_)
    {
#       ifdef FULLDEBUG
        if (debug)
        {
            Info<< "HashTable is empty\n";
        }
#       endif

        return HashTable<T, Key, Hash>::endConstIter_;
    }
    else
    {
        return const_iterator(*this, table_[i], i);
    }
}


template<class T, class Key, class Hash>
inline const typename Foam::HashTable<T, Key, Hash>::const_iterator&
Foam::HashTable<T, Key, Hash>::cend() const
{
    return HashTable<T, Key, Hash>::endConstIter_;
}


template<class T, class Key, class Hash>
inline typename Foam::HashTable<T, Key, Hash>::const_iterator
Foam::HashTable<T, Key, Hash>::begin() const
{
    return this->cbegin();
}


template<class T, class Key, class Hash>
inline const typename Foam::HashTable<T, Key, Hash>::const_iterator&
Foam::HashTable<T, Key, Hash>::end() const
{
    return HashTable<T, Key, Hash>::endConstIter_;
}


// ************************************************************************* //
