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

#include "Switch.H"
#include "error.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// NB: values chosen such that bitwise '&' 0x1 yields the bool value
// INVALID is also evaluates to false, but don't rely on that
const char* Foam::Switch::names[Foam::Switch::INVALID+1] =
{
    "false", "true",
    "off",   "on",
    "no",    "yes",
    "n",     "y",
    "none",  "true",  // is there a reasonable counterpart to "none"?
    "invalid"
};


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

Foam::Switch::switchType Foam::Switch::asEnum(const bool b)
{
    return b ? Switch::TRUE : Switch::FALSE;
}


Foam::Switch::switchType Foam::Switch::asEnum
(
    const std::string& str,
    const bool allowInvalid
)
{
    for (int sw = 0; sw < Switch::INVALID; sw++)
    {
        if (str == names[sw])
        {
            // convert n/y to no/yes (perhaps should deprecate y/n)
            if (sw == Switch::NO_1 || sw == Switch::NONE)
            {
                return Switch::NO;
            }
            else if (sw == Switch::YES_1)
            {
                return Switch::YES;
            }
            else
            {
                return switchType(sw);
            }
        }
    }

    if (!allowInvalid)
    {
        FatalErrorIn("Switch::asEnum(const std::string&)")
            << "unknown switch word " << str << nl
            << abort(FatalError);
    }

    return INVALID;
}


bool Foam::Switch::asBool(const switchType sw)
{
    // relies on (INVALID & 0x1) evaluating to false
    return (sw & 0x1);
}


bool Foam::Switch::asBool
(
    const std::string& str,
    const bool allowInvalid
)
{
    // allow invalid values, but catch after for correct error message
    switchType sw = asEnum(str, true);

    if (sw == Switch::INVALID)
    {
        if (!allowInvalid)
        {
            FatalErrorIn("Switch::asBool(const std::string&)")
                << "unknown switch word " << str << nl
                << abort(FatalError);
        }

        return false;
    }


    return (sw & 0x1);
}


const char* Foam::Switch::asText(const bool b)
{
    return b ? names[Switch::TRUE] : names[Switch::FALSE];
}


const char* Foam::Switch::asText(const switchType sw)
{
    return names[sw];
}


Foam::Switch Foam::Switch::lookupOrAddToDict
(
    const word& name,
    dictionary& dict,
    const Switch& defaultValue
)
{
    return dict.lookupOrAddDefault<Switch>(name, defaultValue);
}


// * * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * //

bool Foam::Switch::readIfPresent(const word& name, const dictionary& dict)
{
    return dict.readIfPresent<Switch>(name, *this);
}


// ************************************************************************* //
