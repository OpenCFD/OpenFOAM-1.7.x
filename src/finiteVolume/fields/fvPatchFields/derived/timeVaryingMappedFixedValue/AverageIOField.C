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

#include "AverageIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
AverageIOField<Type>::AverageIOField
(
    const IOobject& io
)
:
    regIOobject(io),
    pTraits<Type>(readStream(typeName)),
    Field<Type>(readStream(typeName))
{
    close();
}


template<class Type>
AverageIOField<Type>::AverageIOField
(
    const IOobject& io,
    const label size
)
:
    regIOobject(io),
    pTraits<Type>(pTraits<Type>::zero),
    Field<Type>(size)
{}


template<class Type>
AverageIOField<Type>::AverageIOField
(
    const IOobject& io,
    const Type& average,
    const Field<Type>& f
)
:
    regIOobject(io),
    pTraits<Type>(average),
    Field<Type>(f)
{
    if (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    {
        readStream(typeName)
            >> static_cast<Type&>(*this)
            >> static_cast<Field<Type>&>(*this);
        close();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
