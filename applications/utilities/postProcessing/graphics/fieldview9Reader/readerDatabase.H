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
    Foam::readerDatabase

Description
    Singleton caching Foam database and mesh and various. Used in Fv reader
    to keep track of data inbetween callbacks.

SourceFiles
    readerDatabase.C

\*---------------------------------------------------------------------------*/

#ifndef readerDatabase_H
#define readerDatabase_H

#include "wordList.H"
#include "Time.H"
#include "polyMesh.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class fvMeshSubset;
class Time;
class fileName;
class instant;

/*---------------------------------------------------------------------------*\
                           Class readerDatabase Declaration
\*---------------------------------------------------------------------------*/

class readerDatabase
{
    // Private data

        //- Names for protected Fieldview keywords. Gets set at construction
        //  time.
        //  Note: Should be static but this gives problem with construction
        //  order since *this is static as well.
        HashTable<word> fieldviewNames_;

    // Private data

        //- Cached database
        Time* runTimePtr_;

        //- Cached mesh, guaranteed uptodate with runTime.
        fvMeshSubset* meshPtr_;

        //- Empty string or name of current set.
        word setName_;

        //- Cell labels of polyHedra. Uptodate with meshPtr.
        labelList polys_;

        //- All volScalarFields in all time directories
        wordList volScalarNames_;

        //- All volVectorFields ,,
        wordList volVectorNames_;


    // Private Member Functions

        //- Gets cell numbers of all polyHedra
        void getPolyHedra();

        //- Disallow default bitwise copy construct
        readerDatabase(const readerDatabase&);

        //- Disallow default bitwise assignment
        void operator=(const readerDatabase&);


public:

    // Static

        //- Debug flag. Note: uses envvar instead of controlDict since
        //  *this is static as well. Might be initialized before controlDict
        //  read.
        const static bool debug_;


    // Constructors

        //- Construct null
        readerDatabase();


    // Destructor

        ~readerDatabase();


    // Member Functions

        // Access

            const Time& runTime() const;

            const fvMesh& mesh() const;

            const labelList& polys() const;

            const wordList& volScalarNames() const;

            const wordList& volVectorNames() const;

            //- Get fieldview compatible name.
            const word& getFvName(const word& foamName) const;

        // Edit

            //- Create database (if nessecary).
            //  Returns true if new Time created, false if old one reused.
            //  Optional fvMeshSubset using setName.
            bool setRunTime
            (
                const fileName& rootDir,
                const fileName& caseName,
                const word& setName
            );

            //- Forcibly load mesh.
            void loadMesh();

            //- Set time (use this instead of database::setTime), updates
            //  mesh as well and returns mesh update status
            polyMesh::readUpdateState setTime(const instant&, const label);

            //- Set volScalarNames, volVectorNames.
            void setFieldNames(const wordList&, const wordList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
