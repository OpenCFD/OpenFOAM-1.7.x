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

Class
    Foam::keyType

Description
    A class for handling keywords in dictionaries.

    A keyType is the keyword of a dictionary.
    It differs from word in that it accepts patterns (regular expressions).

SourceFiles
    keyType.C
    keyTypeIO.C

\*---------------------------------------------------------------------------*/

#ifndef keyType_H
#define keyType_H

#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;


/*---------------------------------------------------------------------------*\
                          Class keyType Declaration
\*---------------------------------------------------------------------------*/

class keyType
:
    public word
{
    // Private member data

        //- Is the keyType a pattern (regular expression)
        bool isPattern_;

    // Private Member Functions

        //- Disallow assignments where we cannot determine string/word type
        void operator=(const std::string&);

public:


    // Constructors

        //- Construct null
        inline keyType();

        //- Construct as copy
        inline keyType(const keyType&);

        //- Construct as copy of word
        inline keyType(const word&);

        //- Construct as copy of string. Expect it to be regular expression.
        inline keyType(const string&);

        //- Construct as copy of character array
        inline keyType(const char*);

        //- Construct as copy of std::string
        inline keyType(const std::string&, const bool isPattern);

        //- Construct from Istream
        keyType(Istream&);


    // Member functions

        //- Should be treated as a match rather than a literal string
        inline bool isPattern() const;

    // Member operators

        // Assignment

            inline const keyType& operator=(const keyType&);
            inline const keyType& operator=(const word&);

            //- Assign from regular expression.
            inline const keyType& operator=(const string&);
            inline const keyType& operator=(const char*);


    // IOstream operators

        friend Istream& operator>>(Istream&, keyType&);
        friend Ostream& operator<<(Ostream&, const keyType&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "keyTypeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
