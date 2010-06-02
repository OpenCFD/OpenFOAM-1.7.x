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
    Foam::OFstream

Description
    Output to file stream.

SourceFiles
    OFstream.C

\*---------------------------------------------------------------------------*/

#ifndef OFstream_H
#define OFstream_H

#include "OSstream.H"
#include "fileName.H"
#include "className.H"

#include <fstream>
using std::ofstream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OFstream;

/*---------------------------------------------------------------------------*\
                       Class OFstreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- A std::ostream with ability to handle compressed files
class OFstreamAllocator
{
    friend class OFstream;

    ostream* ofPtr_;

    // Constructors

        //- Construct from pathname
        OFstreamAllocator
        (
            const fileName& pathname,
            IOstream::compressionType compression=IOstream::UNCOMPRESSED
        );

    // Destructor

        ~OFstreamAllocator();


public:

    // Member functions

        //- Access to underlying std::ostream
        ostream& stdStream();
};


/*---------------------------------------------------------------------------*\
                           Class OFstream Declaration
\*---------------------------------------------------------------------------*/

class OFstream
:
    private OFstreamAllocator,
    public OSstream
{
    // Private data

        fileName pathname_;


public:

    // Declare name of the class and its debug switch
    ClassName("OFstream");


    // Constructors

        //- Construct from pathname
        OFstream
        (
            const fileName& pathname,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED
        );


    // Destructor

        ~OFstream();


    // Member functions

        // Access

            //- Return the name of the stream
            const fileName& name() const
            {
                return pathname_;
            }

            //- Return non-const access to the name of the stream
            fileName& name()
            {
                return pathname_;
            }


        // Print

            //- Print description of IOstream to Ostream
            void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global predefined null output stream

extern OFstream Snull;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
