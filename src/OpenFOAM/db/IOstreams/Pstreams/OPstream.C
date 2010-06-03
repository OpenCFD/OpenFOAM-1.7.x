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

Description
    Write primitive and binary block from OPstream

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "OPstream.H"
#include "int.H"
#include "token.H"

#include <cctype>

// * * * * * * * * * * * * * Private member functions  * * * * * * * * * * * //

template<class T>
inline void Foam::OPstream::writeToBuffer(const T& t)
{
    writeToBuffer(&t, sizeof(T), sizeof(T));
}


inline void Foam::OPstream::writeToBuffer(const char& c)
{
    if (size_t(buf_.size()) < bufPosition_ + 1U)
    {
        enlargeBuffer(1);
    }

    buf_[bufPosition_] = c;
    bufPosition_ ++;
}


inline void Foam::OPstream::writeToBuffer
(
    const void* data,
    size_t count,
    size_t align
)
{
    label oldPos = bufPosition_;

    if (align > 1)
    {
        // Align bufPosition. Pads bufPosition_ - oldPos characters.
        bufPosition_ = align + ((bufPosition_ - 1) & ~(align - 1));
    }

    if (size_t(buf_.size()) < bufPosition_ + count)
    {
        enlargeBuffer(bufPosition_ - oldPos + count);
    }

    register char* bufPtr = &buf_[bufPosition_];
    register const char* dataPtr = reinterpret_cast<const char*>(data);
    register size_t i = count;
    while (i--) *bufPtr++ = *dataPtr++;

    bufPosition_ += count;
}



// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::OPstream::OPstream
(
    const commsTypes commsType,
    const int toProcNo,
    const label bufSize,
    streamFormat format,
    versionNumber version
)
:
    Pstream(commsType, bufSize),
    Ostream(format, version),
    toProcNo_(toProcNo)
{
    setOpened();
    setGood();

    if (!bufSize)
    {
        buf_.setSize(1000);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::OPstream::write(const token&)
{
    notImplemented("Ostream& OPstream::write(const token&)");
    setBad();
    return *this;
}


Foam::Ostream& Foam::OPstream::write(const char c)
{
    if (!isspace(c))
    {
        writeToBuffer(c);
    }

    return *this;
}


Foam::Ostream& Foam::OPstream::write(const char* str)
{
    word nonWhiteChars(string::validate<word>(str));

    if (nonWhiteChars.size() == 1)
    {
        return write(nonWhiteChars.c_str()[1]);
    }
    else if (nonWhiteChars.size())
    {
        return write(nonWhiteChars);
    }
    else
    {
        return *this;
    }
}


Foam::Ostream& Foam::OPstream::write(const word& str)
{
    write(char(token::WORD));

    size_t len = str.size();
    writeToBuffer(len);
    writeToBuffer(str.c_str(), len + 1, 1);

    return *this;
}


Foam::Ostream& Foam::OPstream::write(const string& str)
{
    write(char(token::STRING));

    size_t len = str.size();
    writeToBuffer(len);
    writeToBuffer(str.c_str(), len + 1, 1);

    return *this;
}


Foam::Ostream& Foam::OPstream::writeQuoted(const std::string& str, const bool)
{
    write(char(token::STRING));

    size_t len = str.size();
    writeToBuffer(len);
    writeToBuffer(str.c_str(), len + 1, 1);

    return *this;
}


Foam::Ostream& Foam::OPstream::write(const label val)
{
    write(char(token::LABEL));
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::OPstream::write(const floatScalar val)
{
    write(char(token::FLOAT_SCALAR));
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::OPstream::write(const doubleScalar val)
{
    write(char(token::DOUBLE_SCALAR));
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::OPstream::write(const char* data, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalErrorIn("Ostream::write(const char*, std::streamsize)")
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    writeToBuffer(data, count, 8);

    return *this;
}


// ************************************************************************* //
