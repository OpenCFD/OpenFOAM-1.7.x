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
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Switch::Switch(Istream& is)
{
    is >> *this;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, Switch& s)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        s.switch_ = Switch::asEnum(bool(t.labelToken()));
    }
    else if (t.isWord())
    {
        // allow invalid values, but catch after for correct error message
        Switch::switchType sw = Switch::asEnum(t.wordToken(), true);

        if (sw == Switch::INVALID)
        {
            is.setBad();
            FatalIOErrorIn("operator>>(Istream&, Switch&)", is)
                << "expected 'true/false', 'on/off' ... found " << t.wordToken()
                << exit(FatalIOError);

            return is;
        }
        else
        {
            s.switch_ = sw;
        }
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, bool/Switch&)", is)
            << "wrong token type - expected bool found " << t
            << exit(FatalIOError);

        return is;
    }


    // Check state of Istream
    is.check("Istream& operator>>(Istream&, Switch&)");

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const Switch& s)
{
    os << Switch::names[s.switch_];
    os.check("Ostream& operator<<(Ostream&, const Switch&)");
    return os;
}


// ************************************************************************* //
