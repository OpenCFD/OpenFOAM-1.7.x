/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
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
    Foam::partialWrite

Description
    Allows some fields/registered objects to be written more often than others.

    Works in the opposite way: deletes at intermediate times all
    but selected fields.

SourceFiles
    partialWrite.C
    IOpartialWrite.H

\*---------------------------------------------------------------------------*/

#ifndef partialWrite_H
#define partialWrite_H

#include "pointFieldFwd.H"
#include "HashSet.H"
#include "DynamicList.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                   Class partialWrite Declaration
\*---------------------------------------------------------------------------*/

class partialWrite
{
protected:

    // Private data

        //- Name of this set of partialWrite
        word name_;

        const objectRegistry& obr_;


        // Read from dictionary

            //- Names of objects to dump always
            HashSet<word> objectNames_;

            //- Write interval for restart dump
            label writeInterval_;


        //- Current dump instance. If reaches writeInterval do a full write.
        label writeInstance_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        partialWrite(const partialWrite&);

        //- Disallow default bitwise assignment
        void operator=(const partialWrite&);


public:

    //- Runtime type information
    TypeName("partialWrite");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        partialWrite
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    //- Destructor
    virtual ~partialWrite();


    // Member Functions

        //- Return name of the partialWrite
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the partialWrite data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Write the partialWrite
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
