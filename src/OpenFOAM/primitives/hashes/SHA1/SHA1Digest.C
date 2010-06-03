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

#include "SHA1Digest.H"
#include "IOstreams.H"

#include <cstring>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

//! @cond fileScope
const char hexChars[] = "0123456789abcdef";
//! @endcond fileScope


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SHA1Digest::SHA1Digest()
{
    clear();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::SHA1Digest::clear()
{
    memset(v_, 0, length);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

bool Foam::SHA1Digest::operator==(const SHA1Digest& rhs) const
{
    for (unsigned i = 0; i < length; ++i)
    {
        if (v_[i] != rhs.v_[i])
        {
            return false;
        }
    }

    return true;
}


bool Foam::SHA1Digest::operator!=(const SHA1Digest& rhs) const
{
    return !operator==(rhs);
}


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const SHA1Digest& dig)
{
    const unsigned char *v = dig.v_;

    for (unsigned i = 0; i < dig.length; ++i)
    {
        os.write(hexChars[((v[i] >> 4) & 0xF)]);
        os.write(hexChars[(v[i] & 0xF)]);
    }

    os.check("Ostream& operator<<(Ostream&, const SHA1Digest&)");
    return os;
}


// ************************************************************************* //
