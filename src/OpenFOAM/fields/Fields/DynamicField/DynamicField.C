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

#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

template<class Type>
const char* const DynamicField<Type>::typeName("DynamicField");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
DynamicField<Type>::DynamicField(Istream& is)
:
    Field<Type>(is),
    capacity_(Field<Type>::size())
{}


template<class Type>
tmp<DynamicField<Type> > DynamicField<Type>::clone() const
{
    return tmp<DynamicField<Type> >(new DynamicField<Type>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void DynamicField<Type>::setSize(const label nElem)
{
    // allocate more capacity?
    if (nElem > capacity_)
    {
        capacity_ = max(nElem, label(1 + capacity_*2));

        Field<Type>::setSize(capacity_);
    }

    // adjust addressed size
    Field<Type>::size(nElem);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * IOstream Operator * * * * * * * * * * * * * //

template<class Type>
Ostream& operator<<(Ostream& os, const DynamicField<Type>& f)
{
    os << static_cast<const Field<Type>&>(f);
    return os;
}


template<class Type>
Ostream& operator<<(Ostream& os, const tmp<DynamicField<Type> >& tf)
{
    os << tf();
    tf.clear();
    return os;
}


template<class Type>
Istream& operator>>(Istream& is, DynamicField<Type>& lst)
{
    is >> static_cast<Field<Type>&>(lst);
    lst.capacity_ = lst.Field<Type>::size();

    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
