/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::fieldMinMax

Description
    Calculates scalar minimim and maximum field values.

    For variables with rank > 0, computes the magnitude of the min/max
    values.

    Data written to the file \<timeDir\>/fieldMinMax.dat


SourceFiles
    fieldMinMax.C
    IOfieldMinMax.H

\*---------------------------------------------------------------------------*/

#ifndef fieldMinMax_H
#define fieldMinMax_H

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "HashSet.H"
#include "OFstream.H"
#include "Switch.H"
#include "pointFieldFwd.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                       Class fieldMinMax Declaration
\*---------------------------------------------------------------------------*/

class fieldMinMax
{
public:

    enum modeType
    {
        mdMag,
        mdCmpt
    };


protected:

    // Protected data

        //- Name of this set of field min/max.
        //  Also used as the name of the output directory.
        word name_;

        const objectRegistry& obr_;

        //- on/off switch
        bool active_;

        //- Switch to send output to Info as well as to file
        Switch log_;

        //- Mode type names
        static const NamedEnum<modeType, 2> modeTypeNames_;

        //- Mode for min/max - only applicable for ranks > 0
        modeType mode_;

        //- Fields to assess min/max
        wordList fieldSet_;


        //- Min/max file ptr
        autoPtr<OFstream> fieldMinMaxFilePtr_;


    // Private Member Functions

        //- If the output file has not been created create it
        void makeFile();

        //- Disallow default bitwise copy construct
        fieldMinMax(const fieldMinMax&);

        //- Disallow default bitwise assignment
        void operator=(const fieldMinMax&);

        //- Output file header information
        virtual void writeFileHeader();


public:

    //- Runtime type information
    TypeName("fieldMinMax");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        fieldMinMax
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false
        );


    // Destructor

        virtual ~fieldMinMax();


    // Member Functions

        //- Return name of the set of field min/max
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the field min/max data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Calculate the field min/max
        template<class Type>
        void calcMinMaxFields(const word& fieldName);

        //- Write the fieldMinMax
        virtual void write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const pointField&)
        {}
};


// Template specialisation for scalar fields
template<>
void fieldMinMax::calcMinMaxFields<scalar>(const word& fieldName);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fieldMinMaxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
