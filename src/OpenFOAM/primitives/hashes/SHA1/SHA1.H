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
    Foam::SHA1

Description
    Functions to compute SHA1 message digest according to the NIST
    specification FIPS-180-1.

    Adapted from the gnulib implementation.

SeeAlso
    Foam::SHA1Digest

SourceFiles
    SHA1I.H
    SHA1.C

\*---------------------------------------------------------------------------*/

#ifndef SHA1_H
#define SHA1_H

#include <string>
#include <cstddef>
#include <stdint.h>    // C++0x uses <cstdint>

#include "SHA1Digest.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Ostream;

// Forward declaration of friend functions and operators
class SHA1;
class SHA1Digest;
Ostream& operator<<(Ostream&, const SHA1&);


/*---------------------------------------------------------------------------*\
                            Class SHA1 Declaration
\*---------------------------------------------------------------------------*/

class SHA1
{
    // Private data

        //- Track if the hashsum has been finalized (added count, etc)
        bool finalized_;

        //- The hash sums
        uint32_t hashsumA_;
        uint32_t hashsumB_;
        uint32_t hashsumC_;
        uint32_t hashsumD_;
        uint32_t hashsumE_;

        //- The total number processed, saved as 64-bit
        uint32_t bufTotal_[2];

        //- The number of elements pending in the buffer
        uint32_t bufLen_;

        //- The input processing buffer
        uint32_t buffer_[32];

    // Private Member Functions

        //- Swap bytes from internal to network (big-endian) order
        static inline uint32_t swapBytes(uint32_t);

        //- Copy the 4-byte value into the memory location pointed to by *dst.
        //  If the architecture allows unaligned access this is equivalent to
        //  *(uint32_t *) cp = val
        static inline void set_uint32(unsigned char *cp, uint32_t);

        //- Process data block-wise, LEN must be a multiple of 64!
        void processBlock(const void *data, size_t len);

        //- Process for the next LEN bytes, LEN need not be a multiple of 64.
        void processBytes(const void *data, size_t len);

        //- Calculate current digest from appended data.
        void calcDigest(SHA1Digest&) const;

public:

    // Constructors

        //- Construct null
        inline SHA1();

        //- Construct and append initial std::string
        explicit inline SHA1(const std::string&);

        //- Construct and append initial string
        explicit inline SHA1(const char*);

    // Member Functions

        //- Reset the hashed data before appending more
        void clear();

        //- Append data for processing
        inline SHA1& append(const char* data, size_t len);

        //- Append string for processing
        inline SHA1& append(const std::string&);

        //- Append string for processing
        inline SHA1& append(const char* str);

        //- Finalized the calculations (normally not needed directly).
        //  Returns false if no bytes were passed for processing
        bool finalize();

        //- Calculate current digest from appended data.
        SHA1Digest digest() const;


    // Member Operators

        //- Equality operator
        inline bool operator==(const SHA1Digest&) const;

        //- Inequality operator
        inline bool operator!=(const SHA1Digest&) const;

        //- Equality operator
        inline bool operator==(const SHA1&) const;

        //- Inequality operator
        inline bool operator!=(const SHA1&) const;

        //- Convert to a digest, calculate current digest from appended data.
        inline operator SHA1Digest() const;

    // Friend Functions

    // Friend Operators

        inline friend Ostream& operator<<(Ostream&, const SHA1&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SHA1I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
