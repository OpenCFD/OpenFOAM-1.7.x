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
    Istream constructor and IOstream operators for word.

\*---------------------------------------------------------------------------*/

#include "word.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word::word(Istream& is)
:
    string()
{
    is >> *this;
}


Foam::Istream& Foam::operator>>(Istream& is, word& w)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isWord())
    {
        w = t.wordToken();
    }
    else if (t.isString())
    {
        // try a bit harder and convert string to word
        w = t.stringToken();
        string::stripInvalid<word>(w);

        // flag empty strings and bad chars as an error
        if (w.empty() || w.size() != t.stringToken().size())
        {
            is.setBad();
            FatalIOErrorIn("operator>>(Istream&, word&)", is)
                << "wrong token type - expected word found non-word characters "
                << t.info()
                << exit(FatalIOError);
            return is;
        }
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, word&)", is)
            << "wrong token type - expected word found "
            << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of IOstream
    is.check("Istream& operator>>(Istream&, word&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const word& w)
{
    os.write(w);
    os.check("Ostream& operator<<(Ostream&, const word&)");
    return os;
}


// ************************************************************************* //
