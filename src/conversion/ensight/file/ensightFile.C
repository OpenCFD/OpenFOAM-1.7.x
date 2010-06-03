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

\*---------------------------------------------------------------------------*/

#include "ensightFile.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::ensightFile::allowUndef_ = false;

Foam::scalar Foam::ensightFile::undefValue_ = Foam::floatScalarVGREAT;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightFile::ensightFile
(
    const fileName& pathname,
    IOstream::streamFormat format
)
:
    OFstream(pathname, format)
{
    // ascii formatting specs
    setf
    (
        ios_base::scientific,
        ios_base::floatfield
    );
    precision(5);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightFile::~ensightFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ensightFile::allowUndef()
{
    return allowUndef_;
}


bool Foam::ensightFile::allowUndef(bool value)
{
    bool old = allowUndef_;
    allowUndef_ = value;
    return old;
}


Foam::scalar Foam::ensightFile::undefValue(const scalar& value)
{
    // enable its use too
    allowUndef_ = true;

    scalar old = undefValue_;
    undefValue_ = value;
    return old;
}


Foam::Ostream& Foam::ensightFile::write
(
    const char* buf,
    std::streamsize count
)
{
    stream().write(buf, count);
    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const string& value)
{
    char buf[80];

    for (string::size_type i = 0; i < 80; ++i)
    {
        buf[i] = 0;
    }

    string::size_type n = value.size();
    if (n >= 80)
    {
        n = 79;
    }

    for (string::size_type i = 0; i < n; ++i)
    {
        buf[i] = value[i];
    }

    if (format() == IOstream::BINARY)
    {
        write
        (
            reinterpret_cast<char const *>(buf),
            sizeof(buf)
        );
    }
    else
    {
        stream() << buf;
    }

    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const label& value)
{
    if (format() == IOstream::BINARY)
    {
        unsigned int ivalue(value);

        write
        (
            reinterpret_cast<char const *>(&ivalue),
            sizeof(ivalue)
        );
    }
    else
    {
        stream().width(10);
        stream() << value;
    }

    return *this;
}


Foam::Ostream& Foam::ensightFile::write
(
    const label& value,
    const label fieldWidth
)
{
    if (format() == IOstream::BINARY)
    {
        unsigned int ivalue(value);

        write
        (
            reinterpret_cast<char const *>(&ivalue),
            sizeof(ivalue)
        );
    }
    else
    {
        stream().width(fieldWidth);
        stream() << value;
    }

    return *this;
}


Foam::Ostream& Foam::ensightFile::write(const scalar& value)
{
    if (format() == IOstream::BINARY)
    {
        float fvalue(value);

        write
        (
            reinterpret_cast<char const *>(&fvalue),
            sizeof(fvalue)
        );
    }
    else
    {
        stream().width(12);
        stream() << value;
    }

    return *this;
}


void Foam::ensightFile::newline()
{
    if (format() == IOstream::ASCII)
    {
        stream() << nl;
    }
}


Foam::Ostream& Foam::ensightFile::writeUndef()
{
    write(undefValue_);
    return *this;
}


Foam::Ostream& Foam::ensightFile::writeKeyword(const string& key)
{
    if (allowUndef_)
    {
        write(key + " undef");
        newline();
        write(undefValue_);
        newline();
    }
    else
    {
        write(key);
        newline();
    }
    return *this;
}


Foam::Ostream& Foam::ensightFile::writeBinaryHeader()
{
    if (format() == IOstream::BINARY)
    {
        write("C Binary");
    }

    return *this;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::string Foam::ensightFile::mask()
{
    char buf[16] = "********";
    return buf;
}


Foam::string Foam::ensightFile::subDir(const label n)
{
    char buf[16];

    sprintf(buf, "%08d", n);
    return buf;
}


// ************************************************************************* //
