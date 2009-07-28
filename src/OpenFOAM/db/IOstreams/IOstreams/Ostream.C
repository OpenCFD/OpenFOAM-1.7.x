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

#include "word.H"
#include "Ostream.H"
#include "token.H"
#include "keyType.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Decrement the indent level
void Foam::Ostream::decrIndent()
{
    if (indentLevel_ == 0)
    {
        cerr<< "Ostream::decrIndent() : attempt to decrement 0 indent level"
            << std::endl;
    }
    else
    {
        indentLevel_--;
    }
}


// Write keyType
// write regular expression as quoted string
// write plain word as word (unquoted)
Foam::Ostream& Foam::Ostream::write(const keyType& kw)
{
    return writeQuoted(kw, kw.isPattern());
}


// Write the keyword followed by appropriate indentation
Foam::Ostream& Foam::Ostream::writeKeyword(const keyType& kw)
{
    indent();
    write(kw);

    label nSpaces = entryIndentation_ - label(kw.size());

    // pattern is surrounded by quotes
    if (kw.isPattern())
    {
        nSpaces -= 2;
    }

    // could also increment by indentSize_ ...
    if (nSpaces < 1)
    {
        nSpaces = 1;
    }

    while (nSpaces--)
    {
        write(char(token::SPACE));
    }

    return *this;
}


// ************************************************************************* //
