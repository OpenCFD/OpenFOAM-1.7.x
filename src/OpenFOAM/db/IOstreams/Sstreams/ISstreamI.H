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

#include "ISstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::ISstream::ISstream
(
    istream& is,
    const string& name,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    Istream(format, version, compression),
    name_(name),
    is_(is)
{
    if (is_.good())
    {
        setOpened();
        setGood();
    }
    else
    {
        setState(is_.rdstate());
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::ISstream& Foam::ISstream::get(char& c)
{
    is_.get(c);
    setState(is_.rdstate());

    if (c == '\n')
    {
        lineNumber_++;
    }

    return *this;
}


inline Foam::ISstream& Foam::ISstream::getLine(string& s)
{
    getline(is_, s);
    setState(is_.rdstate());
    lineNumber_++;

    return *this;
}


inline Foam::ISstream& Foam::ISstream::putback(const char& c)
{
    if (c == '\n')
    {
        lineNumber_--;
    }

    if (!is_.putback(c))
    {
        setBad();
    }

    setState(is_.rdstate());

    return *this;
}


// ************************************************************************* //
