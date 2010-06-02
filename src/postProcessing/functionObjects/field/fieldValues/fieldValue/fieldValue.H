/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::fieldValue

Description
    Base class for field value -based function objects.

SourceFiles
    fieldValue.C

\*---------------------------------------------------------------------------*/

#ifndef fieldValue_H
#define fieldValue_H

#include "Switch.H"
#include "pointFieldFwd.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dictionary;
class objectRegistry;
class fvMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class fieldValue Declaration
\*---------------------------------------------------------------------------*/

class fieldValue
{

protected:

    // Protected data

         //- Name of this fieldValue object
        word name_;

        //- Database this class is registered to
        const objectRegistry& obr_;

        //- Active flag
        bool active_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        //- Name of source object
        word sourceName_;

        //- List of field names to operate on
        wordList fields_;

        //- Output field values flag
        Switch valueOutput_;

        //- Output file pointer
        autoPtr<OFstream> outputFilePtr_;


        // Functions to be over-ridden from IOoutputFilter class

            //- Make the output file
            virtual void makeFile();

            //- Write the output file header
            virtual void writeFileHeader() = 0;

            //- Update mesh
            virtual void updateMesh(const mapPolyMesh&);

            //- Move points
            virtual void movePoints(const Field<point>&);


public:

    //- Run-time type information
    TypeName("fieldValue");


    //- Construct from components
    fieldValue
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        const bool loadFromFiles = false
    );


    //- Destructor
    virtual ~fieldValue();


    // Public member functions

        // Access

            //- Return the name of the geometric source
            inline const word& name() const;

            //- Return the reference to the object registry
            inline const objectRegistry& obr() const;

            //- Return the active flag
            inline bool active() const;

            //- Return the switch to send output to Info as well as to file
            inline const Switch& log() const;

            //- Return the source name
            inline const word& sourceName() const;

            //- Return the list of field names
            inline const wordList& fields() const;

            //- Return the output field values flag
            inline const Switch& valueOutput() const;

            //- Helper function to return the reference to the mesh
            inline const fvMesh& mesh() const;


        // Function object functions

            //- Read from dictionary
            virtual void read(const dictionary& dict);

            //- Write to screen/file
            virtual void write();

            //- Execute
            virtual void execute();

            //- Execute the at the final time-loop, currently does nothing
            virtual void end();

            //- Comnbine fields from all processor domains into single field
            template<class Type>
            tmp<Field<Type> > combineFields
            (
                const tmp<Field<Type> >& field
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fieldValueI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fieldValueTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
