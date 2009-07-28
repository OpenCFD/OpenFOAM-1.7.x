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

#include "timeActivatedExplicitSource.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum
<
    Foam::timeActivatedExplicitSource::volumeType,
    2
>::names[] =
{
    "specific",
    "absolute"
};

const Foam::NamedEnum<Foam::timeActivatedExplicitSource::volumeType, 2>
Foam::timeActivatedExplicitSource::volumeTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeActivatedExplicitSource::timeActivatedExplicitSource
(
    const word& name,
    const fvMesh& mesh,
    const dimensionSet& dims
)
:
    IOdictionary
    (
        IOobject
        (
            name + "Properties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    runTime_(mesh.time()),
    name_(name),
    active_(lookup("active")),
    dimensions_(dims),
    volumeType_(volumeTypeNames_.read(lookup("volumeType"))),
    timeStart_(readScalar(lookup("timeStart"))),
    duration_(readScalar(lookup("duration")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::timeActivatedExplicitSource::mesh() const
{
    return mesh_;
}


const Foam::Time& Foam::timeActivatedExplicitSource::runTime() const
{
    return runTime_;
}


const Foam::word& Foam::timeActivatedExplicitSource::name() const
{
    return name_;
}


const Foam::Switch& Foam::timeActivatedExplicitSource::active() const
{
    return active_;
}


const Foam::dimensionSet& Foam::timeActivatedExplicitSource::dimensions() const
{
    return dimensions_;
}


const Foam::timeActivatedExplicitSource::volumeType&
Foam::timeActivatedExplicitSource::volume() const
{
    return volumeType_;
}


Foam::scalar Foam::timeActivatedExplicitSource::timeStart() const
{
    return timeStart_;
}


Foam::scalar Foam::timeActivatedExplicitSource::duration() const
{
    return duration_;
}


bool Foam::timeActivatedExplicitSource::read()
{
    if (regIOobject::read())
    {
        lookup("active") >> active_;
        if (active_)
        {
            volumeType_ = volumeTypeNames_.read(lookup("volumeType"));
            lookup("timeStart") >> duration_;
            lookup("duration") >> duration_;
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}




// ************************************************************************* //
