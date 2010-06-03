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

#include "CompactListList.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::CompactListList<T>::CompactListList(const List<List<T> >& ll)
:
    offsets_(ll.size())
{
    label sumSize = 0;
    forAll(ll, i)
    {
        sumSize += ll[i].size();
        offsets_[i] = sumSize;
    }

    m_.setSize(sumSize);

    label k = 0;
    forAll(ll, i)
    {
        const List<T>& lli = ll[i];

        forAll(lli, j)
        {
            m_[k++] = lli[j];
        }
    }
}


template<class T>
Foam::CompactListList<T>::CompactListList
(
    const UList<label>& rowSizes
)
:
    offsets_(rowSizes.size())
{
    label sumSize = 0;
    forAll(rowSizes, i)
    {
        sumSize += rowSizes[i];
        offsets_[i] = sumSize;
    }

    m_.setSize(sumSize);
}


template<class T>
Foam::CompactListList<T>::CompactListList
(
    const UList<label>& rowSizes,
    const T& t
)
:
    offsets_(rowSizes.size())
{
    label sumSize = 0;
    forAll(rowSizes, i)
    {
        sumSize += rowSizes[i];
        offsets_[i] = sumSize;
    }

    m_.setSize(sumSize, t);
}


template<class T>
Foam::CompactListList<T>::CompactListList
(
    const Xfer<CompactListList<T> >& lst
)
{
    transfer(lst());
}


template<class T>
Foam::CompactListList<T>::CompactListList
(
    CompactListList<T>& lst,
    bool reUse
)
:
    offsets_(lst.offsets_, reUse),
    m_(lst.m_, reUse)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::CompactListList<T>::setSize(const label nRows)
{
    if (nRows == 0)
    {
        clear();
    }
    if (nRows < offsets_.size())
    {
        offsets_.setSize(nRows);
        m_.setSize(offsets_[nRows - 1]);
    }
    else if (nRows > offsets_.size())
    {
        FatalErrorIn("CompactListList<T>::setSize(const label nRows)")
            << "Cannot be used to extend the list from " << offsets_.size()
            << " to " << nRows << nl
            << "    Please use one of the other setSize member functions"
            << abort(FatalError);
    }
}


template<class T>
void Foam::CompactListList<T>::setSize
(
    const label nRows,
    const label nData
)
{
    offsets_.setSize(nRows);
    m_.setSize(nData);
}


template<class T>
void Foam::CompactListList<T>::setSize
(
    const label nRows,
    const label nData,
    const T& t
)
{
    offsets_.setSize(nRows);
    m_.setSize(nData, t);
}


template<class T>
void Foam::CompactListList<T>::setSize(const UList<label>& rowSizes)
{
    offsets_.setSize(rowSizes.size());

    label sumSize = 0;
    forAll(rowSizes, i)
    {
        sumSize += rowSizes[i];
        offsets_[i] = sumSize;
    }

    m_.setSize(sumSize);
}


template<class T>
Foam::labelList Foam::CompactListList<T>::sizes() const
{
    labelList rowSizes(offsets_.size());

    label prevOffset = 0;
    forAll(offsets_, i)
    {
        rowSizes[i] = offsets_[i]-prevOffset;
        prevOffset = offsets_[i];
    }
    return rowSizes;
}


template<class T>
void Foam::CompactListList<T>::clear()
{
    offsets_.clear();
    m_.clear();
}


template<class T>
void Foam::CompactListList<T>::transfer(CompactListList<T>& a)
{
    offsets_.transfer(a.offsets_);
    m_.transfer(a.m_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
Foam::List<Foam::List<T> > Foam::CompactListList<T>::operator()() const
{
    List<List<T> > ll(offsets_.size());

    label offsetPrev = 0;
    forAll(offsets_, i)
    {
        List<T>& lst = ll[i];

        lst.setSize(offsets_[i] - offsetPrev);

        forAll(lst, j)
        {
            lst[j] = m_[offsetPrev + j];
        }

        offsetPrev = offsets_[i];
    }

    return ll;
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "CompactListListIO.C"

// ************************************************************************* //
