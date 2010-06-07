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
    Foam::surfaceWriter

Description

SourceFiles
    surfaceWriter.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceWriter_H
#define surfaceWriter_H

#include "Field.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "pointField.H"
#include "faceList.H"
#include "fileName.H"

#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type> class surfaceWriter;
template<class Type> class nullSurfaceWriter;
template<class Type> class proxySurfaceWriter;

/*---------------------------------------------------------------------------*\
                        Class surfaceWriter Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class surfaceWriter
{
    //- friendship between writer data types
    template<class Type2> friend class surfaceWriter;

public:

    //- Runtime type information
    TypeName("surfaceWriter");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceWriter,
            word,
            (),
            ()
        );


    // Selectors

        //- Return a reference to the selected surfaceWriter
        static autoPtr<surfaceWriter> New(const word& writeFormat);


    // Constructors

        //- Construct null
        surfaceWriter();


    // Destructor

        virtual ~surfaceWriter();


    // Member Functions

        //- Return true if the surface format supports separate files
        virtual bool separateFiles()
        {
            return false;
        }


        //- Writes single surface geometry to file.
        virtual void write
        (
            const fileName& outputDir,      // <root>/<case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const bool verbose = false
        ) const
        {}


        //- Writes single surface to file. Either one value per vertex or
        //  one value per face (detected by values.size()==faces.size())
        virtual void write
        (
            const fileName& outputDir,      // <root>/<case>/surface/TIME
            const fileName& surfaceName,    // name of surface
            const pointField& points,
            const faceList& faces,
            const fileName& fieldName,      // name of field
            const Field<Type>& values,
            const bool verbose = false
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "surfaceWriter.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
