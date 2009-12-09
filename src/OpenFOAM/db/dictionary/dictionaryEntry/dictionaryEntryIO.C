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

Description
    DictionaryEntry constructor from Istream and Ostream output operator.

\*---------------------------------------------------------------------------*/

#include "keyType.H"
#include "dictionaryEntry.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionaryEntry::dictionaryEntry
(
    const dictionary& parentDict,
    Istream& is
)
:
    entry(keyType(is)),
    dictionary(parentDict, is)
{
    is.fatalCheck
    (
        "dictionaryEntry::dictionaryEntry"
        "(const dictionary& parentDict, Istream&)"
    );
}


Foam::dictionaryEntry::dictionaryEntry
(
    const keyType& key,
    const dictionary& parentDict,
    Istream& is
)
:
    entry(key),
    dictionary(key, parentDict, is)
{
    is.fatalCheck
    (
        "dictionaryEntry::dictionaryEntry"
        "(const keyType&, const dictionary& parentDict, Istream&)"
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dictionaryEntry::write(Ostream& os) const
{
    // write keyword with indent but without trailing spaces
    os.indent();
    os.write(keyword());
    dictionary::write(os);
}


// * * * * * * * * * * * * * Ostream operator  * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const dictionaryEntry& de)
{
    de.write(os);
    return os;
}


#if defined (__GNUC__)
template<>
#endif
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<dictionaryEntry>& ip
)
{
    const dictionaryEntry& e = ip.t_;

    os  << "    dictionaryEntry '" << e.keyword() << "'" << endl;

    return os;
}


// ************************************************************************* //
