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

template<class Type>
inline Foam::SubField<Type>::SubField
(
    const SubList<Type>& list
)
:
    SubList<Type>(list)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const UList<Type>& list
)
:
    SubList<Type>(list, list.size())
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const UList<Type>& list,
    const label subSize
)
:
    SubList<Type>(list, subSize)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const UList<Type>& list,
    const label subSize,
    const label startIndex
)
:
    SubList<Type>(list, subSize, startIndex)
{}


template<class Type>
inline Foam::SubField<Type>::SubField
(
    const SubField<Type>& sfield
)
:
    refCount(),
    SubList<Type>(sfield)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline const Foam::SubField<Type>& Foam::SubField<Type>::null()
{
    //Note: two liner only to work around Icc11 checking
    SubField<Type>* nullPtr = reinterpret_cast< SubField<Type>* >(0);
    return *nullPtr;
}


template<class Type>
inline Foam::tmp<Foam::Field<typename Foam::SubField<Type>::cmptType> >
Foam::SubField<Type>::component
(
    const direction d
) const
{
    return (reinterpret_cast<const Field<Type>&>(*this)).component(d);
}


template<class Type>
inline Foam::tmp<Foam::Field<Type> > Foam::SubField<Type>::T() const
{
    return (reinterpret_cast<const Field<Type>&>(*this)).T();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
inline void Foam::SubField<Type>::operator=(const SubField<Type>& rhs)
{
    UList<Type>::operator=(rhs);
}


template<class Type>
inline Foam::SubField<Type>::operator const Foam::Field<Type>&() const
{
    return *reinterpret_cast< const Field<Type>* >(this);
}


// ************************************************************************* //
