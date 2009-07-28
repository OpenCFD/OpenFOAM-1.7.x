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

#include "pointSourceProperties.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointSourceProperties::pointSourceProperties(Istream& is)
:
    name_("unknownPointSourceName"),
    timeStart_(0.0),
    duration_(0.0),
    location_(point::zero),
    fieldData_()
{
    is.check("pointSourceProperties(Istream&)");

    const dictionaryEntry entry(dictionary::null, is);

    name_ = entry.keyword();
    entry.lookup("timeStart") >> timeStart_;
    entry.lookup("duration") >> duration_;
    entry.lookup("location") >> location_;
    entry.lookup("fieldData") >> fieldData_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, pointSourceProperties& psp)
{
    is.check("Istream& operator>>(Istream&, pointSourceProperties&)");

    const dictionaryEntry entry(dictionary::null, is);

    psp.name_ = entry.keyword();
    entry.lookup("timeStart") >> psp.timeStart_;
    entry.lookup("duration") >> psp.duration_;
    entry.lookup("location") >> psp.location_;
    entry.lookup("fieldData") >> psp.fieldData_;

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const pointSourceProperties& psp)
{
    os.check("Ostream& operator<<(Ostream&, const pointSourceProperties&)");

    os  << psp.name_ << nl << token::BEGIN_BLOCK << nl;
    os.writeKeyword("timeStart") << psp.timeStart_ << token::END_STATEMENT << nl;
    os.writeKeyword("duration") << psp.duration_ << token::END_STATEMENT << nl;
    os.writeKeyword("location") << psp.location_ << token::END_STATEMENT << nl;
    os.writeKeyword("fieldData") << psp.fieldData_ << token::END_STATEMENT << nl;
    os  << token::END_BLOCK << nl;

    os.check("Ostream& operator<<(Ostream&, const pointSourceProperties&)");

    return os;
}



// ************************************************************************* //
