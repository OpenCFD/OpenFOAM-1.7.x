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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::CompactListList<T>::CompactListList()
{}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    const label nRows,
    const label nData
)
:
    offsets_(nRows, 0),
    m_(nData)
{}


template<class T>
inline Foam::CompactListList<T>::CompactListList
(
    const label nRows,
    const label nData,
    const T& t
)
:
    offsets_(nRows, 0),
    m_(nData, t)
{}


template<class T>
inline Foam::autoPtr<Foam::CompactListList<T> >
Foam::CompactListList<T>::clone() const
{
    return autoPtr<CompactListList<T> >(new CompactListList<T>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline const Foam::CompactListList<T>& Foam::CompactListList<T>::null()
{
    return *reinterpret_cast< CompactListList<T>* >(0);
}


template<class T>
inline Foam::label Foam::CompactListList<T>::size() const
{
    return offsets_.size();
}


template<class T>
inline bool Foam::CompactListList<T>::empty() const
{
    return offsets_.empty();
}


template<class T>
inline const Foam::List<Foam::label>& Foam::CompactListList<T>::offsets() const
{
    return offsets_;
}


template<class T>
inline Foam::List<Foam::label>& Foam::CompactListList<T>::offsets()
{
    return offsets_;
}


template<class T>
inline const Foam::List<T>& Foam::CompactListList<T>::m() const
{
    return m_;
}


template<class T>
inline Foam::List<T>& Foam::CompactListList<T>::m()
{
    return m_;
}


template<class T>
inline Foam::label Foam::CompactListList<T>::index
(
    const label i,
    const label j
) const
{
    if (i == 0)
    {
        return j;
    }
    else
    {
        return offsets_[i-1] + j;
    }
}


template<class T>
inline Foam::label Foam::CompactListList<T>::whichRow(const label i) const
{
    if (i < 0 || i >= m_.size())
    {
        FatalErrorIn
        (
            "CompactListList<T>::whichRow(const label) const"
        )   << "Index " << i << " outside 0.." << m_.size()
            << abort(FatalError);
    }

    forAll(offsets_, rowI)
    {
        if (i < offsets_[rowI])
        {
            return rowI;
        }
    }

    return -1;
}


template<class T>
inline Foam::label Foam::CompactListList<T>::whichColumn
(
    const label row,
    const label i
) const
{
    return i - index(row, 0);
}


template<class T>
inline Foam::Xfer<Foam::CompactListList<T> > Foam::CompactListList<T>::xfer()
{
    return xferMove(*this);
}


template<class T>
inline void Foam::CompactListList<T>::resize(const label nRows)
{
    this->setSize(nRows);
}


template<class T>
inline void Foam::CompactListList<T>::resize
(
    const label nRows,
    const label nData
)
{
    this->setSize(nRows, nData);
}


template<class T>
inline void Foam::CompactListList<T>::resize
(
    const label nRows,
    const label nData,
    const T& t
)
{
    this->setSize(nRows, nData, t);
}


template<class T>
inline void Foam::CompactListList<T>::resize(const UList<label>& rowSizes)
{
    this->setSize(rowSizes);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline Foam::UList<T> Foam::CompactListList<T>::operator[]
(
    const label i
)
{
    if (i == 0)
    {
        return UList<T>(m_.begin(), offsets_[i]);
    }
    else
    {
        return UList<T>(&m_[offsets_[i-1]], offsets_[i] - offsets_[i-1]);
    }
}


template<class T>
inline const Foam::UList<T> Foam::CompactListList<T>::operator[]
(
    const label i
) const
{
    if (i == 0)
    {
        return UList<T>(m_.begin(), offsets_[i]);
    }
    else
    {
        return UList<T>(&m_[offsets_[i-1]], offsets_[i] - offsets_[i-1]);
    }
}


template<class T>
inline T& Foam::CompactListList<T>::operator()
(
    const label i,
    const label j
)
{
    return m_[index(i, j)];
}


template<class T>
inline const T& Foam::CompactListList<T>::operator()
(
    const label i,
    const label j
) const
{
    return m_[index(i, j)];
}


template<class T>
inline void Foam::CompactListList<T>::operator=(const T& t)
{
    m_ = t;
}


// ************************************************************************* //
