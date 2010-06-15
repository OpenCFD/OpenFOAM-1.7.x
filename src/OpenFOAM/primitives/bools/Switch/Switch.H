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
    Foam::Switch

Description
    A simple wrapper around bool so that it can be read as a word:
    true/false, on/off, yes/no or y/n or none.

SourceFiles
    Switch.C
    SwitchIO.C

\*---------------------------------------------------------------------------*/

#ifndef Switch_H
#define Switch_H

#include "bool.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class Switch;

Istream& operator>>(Istream&, Switch&);
Ostream& operator<<(Ostream&, const Switch&);

class dictionary;

/*---------------------------------------------------------------------------*\
                           Class Switch Declaration
\*---------------------------------------------------------------------------*/

class Switch
{
private:

    // Private data

        //- The logic and enumerated text representation stored as a single byte
        unsigned char switch_;

public:

    // Public data types

        //- The various text representations for a switch value.
        //  These also correspond to the entries in names.
        enum switchType
        {
            FALSE = 0,  TRUE  = 1,
            OFF   = 2,  ON    = 3,
            NO    = 4,  YES   = 5,
            NO_1  = 6,  YES_1 = 7,
            NONE  = 8,  PLACEHOLDER = 9,
            INVALID
        };


    // Static data members

        //- The set of names corresponding to the switchType enumeration
        //  Includes an extra entry for "invalid".
        static const char* names[INVALID+1];


    // Static Member Functions

        //- Return a switchType representation of a bool
        static switchType asEnum(const bool);

        //- Return a switchType representation of a word
        //  Optionally allow bad words, and catch the error elsewhere
        static switchType asEnum
        (
            const std::string&,
            const bool allowInvalid=false
        );

        //- Return a bool representation of a switchType
        static bool asBool(const switchType);

        //- Return a bool representation of a word
        //  Optionally allow bad words, and catch the error elsewhere
        static bool asBool
        (
            const std::string&,
            const bool allowInvalid=false
        );

        //- Return a text representation of a bool value
        static const char* asText(const bool);

        //- Return a text representation of a switchType
        static const char* asText(const switchType);


    // Constructors

        //- Construct null as false
        Switch()
        :
            switch_(Switch::FALSE)
        {}

        //- Construct from bool
        Switch(const bool value)
        :
            switch_(asEnum(value))
        {}

        //- Construct from integer values (treats integer as bool value)
        Switch(const int value)
        :
            switch_(asEnum(bool(value)))
        {}

        //- Construct from std::string, string, word
        Switch(const std::string& value)
        :
            switch_(asEnum(value))
        {}

        //- Construct from character array
        Switch(const char* value)
        :
            switch_(asEnum(std::string(value)))
        {}

        //- Construct from Istream
        Switch(Istream& is);

        //- Construct from dictionary, supplying default value so that if the
        //  value is not found, it is added into the dictionary.
        static Switch lookupOrAddToDict
        (
            const word&,
            dictionary&,
            const Switch& defaultValue = false
        );


    // Member Operators

        //- Conversion to bool
        operator bool() const
        {
            return (switch_ & 0x1);
        }

        //- Assignment from bool
        const Switch& operator=(const bool b)
        {
            switch_ = (b ? Switch::TRUE : Switch::FALSE);
            return *this;
        }


    // Member fuctions

        //- Update the value of the Switch if it is found in the dictionary
        bool readIfPresent(const word&, const dictionary&);


    // IOstream Operators

        friend Istream& operator>>(Istream&, Switch&);
        friend Ostream& operator<<(Ostream&, const Switch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
