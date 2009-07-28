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

#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionaryEntry::dictionaryEntry
(
    const keyType& key,
    const dictionary& parentDict,
    const dictionary& dict
)
:
    entry(key),
    dictionary(parentDict, dict)
{}


Foam::dictionaryEntry::dictionaryEntry
(
    const dictionary& parentDict,
    const dictionaryEntry& dictEnt
)
:
    entry(dictEnt),
    dictionary(parentDict, dictEnt)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::dictionaryEntry::startLineNumber() const
{
    if (size())
    {
        return first()->startLineNumber();
    }
    else
    {
        return -1;
    }
}

Foam::label Foam::dictionaryEntry::endLineNumber() const
{
    if (size())
    {
        return last()->endLineNumber();
    }
    else
    {
        return -1;
    }
}


Foam::ITstream& Foam::dictionaryEntry::stream() const
{
    FatalIOErrorIn("ITstream& primitiveEntry::stream() const", *this)
        << "Attempt to return dictionary entry as a primitive"
        << abort(FatalIOError);

    return lookup("");
}


const Foam::dictionary& Foam::dictionaryEntry::dict() const
{
    return *this;
}


Foam::dictionary& Foam::dictionaryEntry::dict()
{
    return *this;
}


// ************************************************************************* //
