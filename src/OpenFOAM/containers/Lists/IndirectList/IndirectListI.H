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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
inline Foam::IndirectList<T>::IndirectList
(
    const UList<T>& completeList,
    const UList<label>& addr
)
:
    completeList_(const_cast<UList<T>&>(completeList)),
    addressing_(addr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::IndirectList<T>::size() const
{
    return addressing_.size();
}


template<class T>
inline bool Foam::IndirectList<T>::empty() const
{
    return addressing_.empty();
}


template<class T>
inline const Foam::UList<T>& Foam::IndirectList<T>::completeList() const
{
    return completeList_;
}


template<class T>
inline const Foam::List<Foam::label>& Foam::IndirectList<T>::addressing() const
{
    return addressing_;
}


template<class T>
inline void Foam::IndirectList<T>::resetAddressing
(
    const UList<label>& addr
)
{
    addressing_ = addr;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline Foam::List<T> Foam::IndirectList<T>::operator()() const
{
    List<T> result(size());

    forAll(*this, i)
    {
        result[i] = operator[](i);
    }

    return result;
}


template<class T>
inline T& Foam::IndirectList<T>::operator[](const label i)
{
    return completeList_[addressing_[i]];
}


template<class T>
inline const T& Foam::IndirectList<T>::operator[](const label i) const
{
    return completeList_[addressing_[i]];
}


template<class T>
inline void Foam::IndirectList<T>::operator=(const UList<T>& ae)
{
    if (addressing_.size() != ae.size())
    {
        FatalErrorIn("IndirectList<T>::operator=(const UList<T>&)")
            << "Addressing and list of addressed elements "
               "have different sizes: "
            << addressing_.size() << " " << ae.size()
            << abort(FatalError);
    }

    forAll(addressing_, i)
    {
        completeList_[addressing_[i]] = ae[i];
    }
}


template<class T>
inline void Foam::IndirectList<T>::operator=(const T& t)
{
    forAll(addressing_, i)
    {
        completeList_[addressing_[i]] = t;
    }
}


// ************************************************************************* //
