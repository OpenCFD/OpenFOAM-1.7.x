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
    Reads an uint from an input stream, for a given version
    number and File format. If an ascii File is being read,
    then the line numbers are counted and an erroneous read
    ised.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "uint.H"
#include "IOstreams.H"

#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const unsigned int val)
{
    std::ostringstream buf;
    buf << val;
    return buf.str();
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, unsigned int& i)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        i = static_cast<unsigned int>(t.labelToken());
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, unsigned int&)", is)
            << "wrong token type - expected unsigned int found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, unsigned int&)");

    return is;
}


unsigned int Foam::readUint(Istream& is)
{
    unsigned int val;
    is >> val;

    return val;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const unsigned int i)
{
    os.write(label(i));
    os.check("Ostream& operator<<(Ostream&, const unsigned int)");
    return os;
}


// ************************************************************************* //
