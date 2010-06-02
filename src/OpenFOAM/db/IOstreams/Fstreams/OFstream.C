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

#include "OFstream.H"
#include "OSspecific.H"
#include "gzstream.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::OFstream, 0);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::OFstreamAllocator::OFstreamAllocator
(
    const fileName& pathname,
    IOstream::compressionType compression
)
:
    ofPtr_(NULL)
{
    if (pathname.empty())
    {
        if (OFstream::debug)
        {
            Info<< "OFstreamAllocator::OFstreamAllocator(const fileName&) : "
                   "cannot open null file " << endl;
        }
    }

    if (compression == IOstream::COMPRESSED)
    {
        // get identically named uncompressed version out of the way
        if (isFile(pathname, false))
        {
            rm(pathname);
        }

        ofPtr_ = new ogzstream((pathname + ".gz").c_str());
    }
    else
    {
        // get identically named compressed version out of the way
        if (isFile(pathname + ".gz", false))
        {
            rm(pathname + ".gz");
        }

        ofPtr_ = new ofstream(pathname.c_str());
    }
}


Foam::OFstreamAllocator::~OFstreamAllocator()
{
    delete ofPtr_;
}


std::ostream& Foam::OFstreamAllocator::stdStream()
{
    if (!ofPtr_)
    {
        FatalErrorIn("OFstreamAllocator::stdStream()")
            << "No stream allocated." << abort(FatalError);
    }
    return *ofPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::OFstream::OFstream
(
    const fileName& pathname,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    OFstreamAllocator(pathname, compression),
    OSstream(*ofPtr_, "OFstream.sinkFile_", format, version, compression),
    pathname_(pathname)
{
    setClosed();
    setState(ofPtr_->rdstate());

    if (!good())
    {
        if (debug)
        {
            Info<< "IFstream::IFstream(const fileName&,"
                   "streamFormat format=ASCII,"
                   "versionNumber version=currentVersion) : "
                   "could not open file for input\n"
                   "in stream " << info() << Foam::endl;
        }

        setBad();
    }
    else
    {
        setOpened();
    }

    lineNumber_ = 1;
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::OFstream::~OFstream()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::OFstream::print(Ostream& os) const
{
    os  << "    OFstream: ";
    OSstream::print(os);
}


// ************************************************************************* //
