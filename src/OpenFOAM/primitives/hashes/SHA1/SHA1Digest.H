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

Class
    Foam::SHA1Digest

Description
    The SHA1 message digest.

SeeAlso
    Foam::SHA1

SourceFiles
    SHA1Digest.C

\*---------------------------------------------------------------------------*/

#ifndef SHA1Digest_H
#define SHA1Digest_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;

// Forward declaration of friend functions and operators
class SHA1;
class SHA1Digest;
Ostream& operator<<(Ostream&, const SHA1Digest&);


/*---------------------------------------------------------------------------*\
                        Class SHA1Digest Declaration
\*---------------------------------------------------------------------------*/

class SHA1Digest
{
public:
    friend class SHA1;

    //- The length of the digest contents
    static const unsigned length = 20;

    //- Construct a zero digest
    SHA1Digest();

    //- Reset the digest to zero
    void clear();

    //- Equality operator
    bool operator==(const SHA1Digest&) const;

    //- Inequality operator
    bool operator!=(const SHA1Digest&) const;

    friend Ostream& operator<<(Ostream&, const SHA1Digest&);

private:

    //- The digest contents
    unsigned char v_[length];
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
