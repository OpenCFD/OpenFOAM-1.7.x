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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointSourceProperties::pointSourceProperties()
:
    name_("unknownPointSourceName"),
    timeStart_(0.0),
    duration_(0.0),
    location_(point::zero),
    fieldData_()
{}


Foam::pointSourceProperties::pointSourceProperties(const dictionary& dict)
:
    name_(dict.name().name()),
    timeStart_(readScalar(dict.lookup("timeStart"))),
    duration_(readScalar(dict.lookup("duration"))),
    location_(dict.lookup("location")),
    fieldData_(dict.lookup("fieldData"))
{}


Foam::pointSourceProperties::pointSourceProperties
(
    const pointSourceProperties& psp
)
:
    name_(psp.name_),
    timeStart_(psp.timeStart_),
    duration_(psp.duration_),
    location_(psp.location_),
    fieldData_(psp.fieldData_)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pointSourceProperties::operator=(const pointSourceProperties& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "pointSourceProperties::operator=(const pointSourceProperties&)"
        )   << "Attempted assignment to self" << nl
            << abort(FatalError);
    }

    // Set updated values
    name_ = rhs.name_;
    timeStart_ = rhs.timeStart_;
    duration_ = rhs.duration_;
    location_ = rhs.location_;
    fieldData_ = rhs.fieldData_;}



// ************************************************************************* //
