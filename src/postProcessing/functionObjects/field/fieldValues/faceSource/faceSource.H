/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    Foam::fieldValues::faceSource

Description
    Face source variant of field value function object. Values of user-
    specified fields reported for collections of faces.

    faceObj1                        // Name also used to identify output folder
    {
        type            faceSource;
        functionObjectLibs ("libfieldValueFunctionObjects.so");
        enabled         true;
        outputControl   outputTime;
        log             true;       // log to screen?
        valueOutput     true;       // Write values at run-time output times?
        source          faceZone;   // Type of face source: faceZone, patch
        sourceName      f0;
        operation       sum;
        fields
        (
            p
            phi
            U
        );
    }

    where operation is one of:
      - none
      - sum
      - areaAverage
      - areaIntegrate
      - weightedAverage

SourceFiles
    faceSource.C

\*---------------------------------------------------------------------------*/

#ifndef faceSource_H
#define faceSource_H

#include "NamedEnum.H"
#include "fieldValue.H"
#include "labelList.H"
#include "surfaceFieldsFwd.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                         Class faceSource Declaration
\*---------------------------------------------------------------------------*/

class faceSource
:
    public fieldValue
{

public:

    // Public data types

        //- Source type enumeration
        enum sourceType
        {
            stFaceZone,
            stPatch
        };

        //- Source type names
        static const NamedEnum<sourceType, 2> sourceTypeNames_;


        //- Operation type enumeration
        enum operationType
        {
            opNone,
            opSum,
            opAreaAverage,
            opAreaIntegrate,
            opWeightedAverage
        };

        //- Operation type names
        static const NamedEnum<operationType, 5> operationTypeNames_;


private:

    // Private member functions

        //- Set faces to evaluate based on a face zone
        void setFaceZoneFaces();

        //- Set faces to evaluate based on a patch
        void setPatchFaces();


protected:

    // Protected data

        //- Source type
        sourceType source_;

        //- Operation to apply to values
        operationType operation_;

        //- Local list of face IDs
        labelList faceId_;

        //- Local list of patch ID per face
        labelList facePatchId_;

        //- List of +1/-1 representing face flip map
        labelList flipMap_;

        //- Weight field name - only used for opWeightedAverage mode
        word weightFieldName_;


    // Protected member functions

        //- Initialise, e.g. face addressing
        void initialise(const dictionary& dict);

        //- Return true if the field name is valid
        template<class Type>
        bool validField(const word& fieldName) const;

        //- Return field values by looking up field name
        template<class Type>
        tmp<Field<Type> > setFieldValues(const word& fieldName) const;

        //- Apply the 'operation' to the values
        template<class Type>
        Type processValues
        (
            const Field<Type>& values,
            const scalarField& magSf,
            const scalarField& weightField
        ) const;

        //- Output file header information
        virtual void writeFileHeader();


public:

    //- Run-time type information
    TypeName("faceSource");


    //- Construct from components
    faceSource
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        const bool loadFromFiles = false
    );


    //- Destructor
    virtual ~faceSource();


    // Public member functions

        // Access

            //- Return the source type
            inline const sourceType& source() const;

            //- Return the local list of face IDs
            inline const labelList& faceId() const;

            //- Return the local list of patch ID per face
            inline const labelList& facePatch() const;

            //- Return the list of +1/-1 representing face flip map
            inline const labelList& flipMap() const;


        // Function object functions

            //- Read from dictionary
            virtual void read(const dictionary&);

            //- Calculate and write
            virtual void write();

            //- Templated helper function to output field values
            template<class Type>
            bool writeValues(const word& fieldName);

            //- Filter a surface field according to faceIds
            template<class Type>
            tmp<Field<Type> > filterField
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>& field
            ) const;

            //- Filter a volume field according to faceIds
            template<class Type>
            tmp<Field<Type> > filterField
            (
                const GeometricField<Type, fvPatchField, volMesh>& field
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "faceSourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faceSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
