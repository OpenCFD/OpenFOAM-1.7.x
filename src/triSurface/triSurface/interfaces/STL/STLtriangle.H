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
    Foam::STLtriangle

Description
    A triangle representation for STL files.

SourceFiles
    STLtriangleI.H

\*---------------------------------------------------------------------------*/

#ifndef STLtriangle_H
#define STLtriangle_H

#include "STLpoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class STLtriangle Declaration
\*---------------------------------------------------------------------------*/

class STLtriangle
{
    // Private data

        STLpoint normal_, a_, b_, c_;
        unsigned short region_;


public:

    // Constructors

        //- Construct null
        inline STLtriangle();

        //- Construct from components
        inline STLtriangle
        (
            const STLpoint& normal,
            const STLpoint& a,
            const STLpoint& b,
            const STLpoint& c,
            unsigned short region
        );

        //- Construct from istream (read binary)
        inline STLtriangle(istream&);


    // Member Functions

        // Access

            inline const STLpoint& a() const;
            inline const STLpoint& b() const;
            inline const STLpoint& c() const;
            inline unsigned short region() const;

        // Read

            //- Read from istream (binary)
            inline void read(istream&);

        // Write

            //- Write to istream (binary)
            inline void write(ostream&);


    // Ostream operator

        inline friend Ostream& operator<<(Ostream&, const STLtriangle&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "STLtriangleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
