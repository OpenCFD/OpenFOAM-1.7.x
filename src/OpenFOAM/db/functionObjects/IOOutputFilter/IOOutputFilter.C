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

\*---------------------------------------------------------------------------*/

#include "IOOutputFilter.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class OutputFilter>
Foam::IOOutputFilter<OutputFilter>::IOOutputFilter
(
    const word& outputFilterName,
    const objectRegistry& obr,
    const fileName& dictName,
    const IOobject::readOption rOpt,
    const bool readFromFiles
)
:
    IOdictionary
    (
        IOobject
        (
            dictName,
            obr.time().system(),
            obr,
            rOpt,
            IOobject::NO_WRITE
        )
    ),
    OutputFilter(outputFilterName, obr, *this, readFromFiles)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class OutputFilter>
Foam::IOOutputFilter<OutputFilter>::~IOOutputFilter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class OutputFilter>
bool Foam::IOOutputFilter<OutputFilter>::read()
{
    if (regIOobject::read())
    {
        OutputFilter::read(*this);
        return true;
    }
    else
    {
        return false;
    }
}


template<class OutputFilter>
void Foam::IOOutputFilter<OutputFilter>::write()
{
    OutputFilter::write();
}


// ************************************************************************* //
