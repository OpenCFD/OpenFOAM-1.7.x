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
    Foam::writeRegisteredObject

Description
    Takes over the writing of registered IO objects

SourceFiles
    writeRegisteredObject.C
    IOwriteRegisteredObject.H

\*---------------------------------------------------------------------------*/

#ifndef writeRegisteredObject_H
#define writeRegisteredObject_H

#include "pointFieldFwd.H"
#include "wordList.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                   Class writeRegisteredObject Declaration
\*---------------------------------------------------------------------------*/

class writeRegisteredObject
{
protected:

    // Private data

        //- Name of this set of writeRegisteredObject
        word name_;

        const objectRegistry& obr_;

        //- On/off switch
        bool active_;


        // Read from dictionary

            //- Names of objects to control
            wordList objectNames_;


    // Private Member Functions


        //- Disallow default bitwise copy construct
        writeRegisteredObject(const writeRegisteredObject&);

        //- Disallow default bitwise assignment
        void operator=(const writeRegisteredObject&);


public:

    //- Runtime type information
    TypeName("writeRegisteredObject");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        writeRegisteredObject
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~writeRegisteredObject();


    // Member Functions

        //- Return name of the writeRegisteredObject
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the writeRegisteredObject data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Write the writeRegisteredObject
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
