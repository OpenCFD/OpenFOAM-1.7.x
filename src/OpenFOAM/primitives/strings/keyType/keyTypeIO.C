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

Description
    Istream constructor and IOstream operators for word.

\*---------------------------------------------------------------------------*/

#include "keyType.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::keyType::keyType(Istream& is)
:
    word()
{
    is >> *this;
}


Foam::Istream& Foam::operator>>(Istream& is, keyType& w)
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
        // Assign from string. Sets regular expression.
        w = t.stringToken();
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, keyType&)", is)
            << "wrong token type - expected word or string found "
            << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of IOstream
    is.check("Istream& operator>>(Istream&, keyType&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const keyType& w)
{
    os.write(w);
    os.check("Ostream& operator<<(Ostream&, const keyType&)");
    return os;
}


// ************************************************************************* //
