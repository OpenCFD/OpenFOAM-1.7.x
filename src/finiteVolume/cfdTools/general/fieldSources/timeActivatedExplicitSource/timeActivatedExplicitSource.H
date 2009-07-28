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

Class
    Foam::timeActivatedExplicitSource

Description
    Base class for field sources. Provides:
    - name
    - references to mesh and time
    - dimensions
    - volume type
    - startTime
    - duration
    - read (from <sourceName>Properties dictionary)

SourceFiles
    timeActivatedExplicitSource.C

\*---------------------------------------------------------------------------*/

#ifndef timeActivatedExplicitSource_H
#define timeActivatedExplicitSource_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "fvMesh.H"
#include "Time.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class timeActivatedExplicitSource Declaration
\*---------------------------------------------------------------------------*/

class timeActivatedExplicitSource
:
    public IOdictionary
{
public:

    enum volumeType
    {
        vtSpecific,
        vtAbsolute
    };

    static const NamedEnum<volumeType, 2> volumeTypeNames_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        timeActivatedExplicitSource(const timeActivatedExplicitSource&);

        //- Disallow default bitwise assignment
        void operator=(const timeActivatedExplicitSource&);


protected:

    // Protected data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to time database
        const Time& runTime_;


        // Source properties

            //- Name of the source
            word name_;

            //- Active flag
            Switch active_;

            //- Dimensions
            const dimensionSet dimensions_;

            //- Volume type
            volumeType volumeType_;

            //- Time start [s]
            scalar timeStart_;

            //- Duration [s]
            scalar duration_;


public:

    // Constructors

        //- Construct from explicit source name and mesh
        timeActivatedExplicitSource
        (
            const word&,
            const fvMesh&,
            const dimensionSet&
        );


    // Member Functions

        // Access

            //- Return the reference to the mesh
            virtual const fvMesh& mesh() const;

            //- Return the reference to the time database
            virtual const Time& runTime() const;

            //- Return the source name
            virtual const word& name() const;

            //- Return the active flag
            virtual const Switch& active() const;

            //- Return the dimensions
            virtual const dimensionSet& dimensions() const;

            //- Return the volume type
            virtual const volumeType& volume() const;

            //- Return the start time
            virtual scalar timeStart() const;

            //- Return the duration
            virtual scalar duration() const;


        //- Read properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
