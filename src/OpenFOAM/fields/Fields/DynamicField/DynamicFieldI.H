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

#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
DynamicField<Type>::DynamicField()
:
    Field<Type>(),
    capacity_(0)
{}


template<class Type>
DynamicField<Type>::DynamicField(const label size)
:
    Field<Type>(size),
    capacity_(Field<Type>::size())
{
    Field<Type>::size(0);
}


template<class Type>
inline Foam::DynamicField<Type>::DynamicField
(
    const UList<Type>& lst
)
:
    Field<Type>(lst),
    capacity_(Field<Type>::size())
{}


template<class Type>
inline Foam::DynamicField<Type>::DynamicField
(
    const Xfer<List<Type> >& lst
)
:
    Field<Type>(lst),
    capacity_(Field<Type>::size())
{}


template<class Type>
DynamicField<Type>::DynamicField
(
    const UList<Type>& mapF,
    const labelList& mapAddressing
)
:
    Field<Type>(mapF, mapAddressing),
    capacity_(Field<Type>::size())
{}


template<class Type>
DynamicField<Type>::DynamicField
(
    const UList<Type>& mapF,
    const labelListList& mapAddressing,
    const scalarListList& weights
)
:
    Field<Type>(mapF, mapAddressing, weights),
    capacity_(Field<Type>::size())
{}


//- Construct by mapping from the given field
template<class Type>
DynamicField<Type>::DynamicField
(
    const UList<Type>& mapF,
    const FieldMapper& map
)
:
    DynamicField<Type>(mapF, map),
    capacity_(Field<Type>::size())
{}


template<class Type>
DynamicField<Type>::DynamicField(const DynamicField<Type>& f)
:
    Field<Type>(f),
    capacity_(Field<Type>::size())
{}


template<class Type>
DynamicField<Type>::DynamicField(DynamicField<Type>& f, bool reUse)
:
    Field<Type>(f, reUse),
    capacity_(Field<Type>::size())
{}


template<class Type>
DynamicField<Type>::DynamicField(const Xfer<DynamicField<Type> >& f)
:
    Field<Type>(f),
    capacity_(Field<Type>::size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::label DynamicField<Type>::capacity() const
{
    return capacity_;
}


template<class Type>
void DynamicField<Type>::append(const Type& t)
{
    label elemI = Field<Type>::size();
    setSize(elemI + 1);

    this->operator[](elemI) = t;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void DynamicField<Type>::operator=(const DynamicField<Type>& rhs)
{
    if (this == &rhs)
    {
        FatalErrorIn("DynamicField<Type>::operator=(const DynamicField<Type>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    Field<Type>::operator=(rhs);
    capacity_ = Field<Type>::size();
}


template<class Type>
void DynamicField<Type>::operator=(const UList<Type>& rhs)
{
    Field<Type>::operator=(rhs);
    capacity_ = Field<Type>::size();
}


template<class Type>
void DynamicField<Type>::operator=(const tmp<DynamicField>& rhs)
{
    if (this == &(rhs()))
    {
        FatalErrorIn("DynamicField<Type>::operator=(const tmp<DynamicField>&)")
            << "attempted assignment to self"
            << abort(FatalError);
    }

    // This is dodgy stuff, don't try it at home.
    DynamicField* fieldPtr = rhs.ptr();
    List<Type>::transfer(*fieldPtr);
    delete fieldPtr;
    capacity_ = Field<Type>::size();
}


template<class Type>
Type& DynamicField<Type>::operator[](const label i)
{
    return Field<Type>::operator[](i);
}


template<class Type>
const Type& DynamicField<Type>::operator[](const label i) const
{
    return Field<Type>::operator[](i);
}


// * * * * * * * * * * * * * * * IOstream Operator * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
