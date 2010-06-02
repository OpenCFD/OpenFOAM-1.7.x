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

#include "Constant.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Constant<Type>::Constant(const word& entryName, Istream& is)
:
    DataEntry<Type>(entryName),
    value_(is)
{}


template<class Type>
Foam::Constant<Type>::Constant(const Constant<Type>& cnst)
:
    DataEntry<Type>(cnst),
    value_(cnst.value_)
{}


template<>
Foam::Constant<Foam::label>::Constant(const word& entryName, Istream& is)
:
    DataEntry<label>(entryName),
    value_(readLabel(is))
{}


template<>
Foam::Constant<Foam::scalar>::Constant(const word& entryName, Istream& is)
:
    DataEntry<scalar>(entryName),
    value_(readScalar(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Constant<Type>::~Constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Constant<Type>::value(const scalar x) const
{
    return value_;
}


template<class Type>
Type Foam::Constant<Type>::integrate(const scalar x1, const scalar x2) const
{
    return (x2 - x1)*value_;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "ConstantIO.C"


// ************************************************************************* //
