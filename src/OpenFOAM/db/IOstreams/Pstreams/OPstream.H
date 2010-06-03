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
    Foam::OPstream

Description
    Output inter-processor communications stream.

SourceFiles
    OPstream.C

\*---------------------------------------------------------------------------*/

#include "Pstream.H"

#ifndef OPstream_H
#define OPstream_H

#include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class OPstream Declaration
\*---------------------------------------------------------------------------*/

class OPstream
:
    public Pstream,
    public Ostream
{
    // Private data

        int toProcNo_;


    // Private member functions

        //- Write a T to the transfer buffer
        template<class T>
        inline void writeToBuffer(const T&);

        //- Write a char to the transfer buffer
        inline void writeToBuffer(const char&);

        //- Write data to the transfer buffer
        inline void writeToBuffer(const void* data, size_t count, size_t align);


public:

    // Constructors

        //- Construct given process index to send to and optional buffer size,
        //  write format and IO version
        OPstream
        (
            const commsTypes commsType,
            const int toProcNo,
            const label bufSize = 0,
            streamFormat format=BINARY,
            versionNumber version=currentVersion
        );


    // Destructor

        ~OPstream();


    // Member functions

        // Inquiry

            //- Return flags of output stream
            ios_base::fmtflags flags() const
            {
                return ios_base::fmtflags(0);
            }


        // Write functions

            //- Write given buffer to given processor
            static bool write
            (
                const commsTypes commsType,
                const int toProcNo,
                const char* buf,
                const std::streamsize bufSize
            );

            //- Non-blocking writes: wait until all have finished.
            static void waitRequests();

            //- Non-blocking writes: has request i finished?
            static bool finishedRequest(const label i);

            //- Write next token to stream
            Ostream& write(const token&);

            //- Write character
            Ostream& write(const char);

            //- Write character string
            Ostream& write(const char*);

            //- Write word
            Ostream& write(const word&);

            //- Write string
            Ostream& write(const string&);

            //- Write std::string surrounded by quotes.
            //  Optional write without quotes.
            Ostream& writeQuoted
            (
                const std::string&,
                const bool quoted=true
            );

            //- Write label
            Ostream& write(const label);

            //- Write floatScalar
            Ostream& write(const floatScalar);

            //- Write doubleScalar
            Ostream& write(const doubleScalar);

            //- Write binary block
            Ostream& write(const char*, std::streamsize);

            //- Add indentation characters
            void indent()
            {}


        // Stream state functions

            //- Flush stream
            void flush()
            {}

            //- Add newline and flush stream
            void endl()
            {}

            //- Get width of output field
            int width() const
            {
                return 0;
            }

            //- Set width of output field (and return old width)
            int width(const int)
            {
                 return 0;
            }

            //- Get precision of output field
            int precision() const
            {
                 return 0;
            }

            //- Set precision of output field (and return old precision)
            int precision(const int)
            {
                 return 0;
            }


        // Edit

            //- Set flags of stream
            ios_base::fmtflags flags(const ios_base::fmtflags)
            {
                return ios_base::fmtflags(0);
            }


        // Print

            //- Print description of IOstream to Ostream
            void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
