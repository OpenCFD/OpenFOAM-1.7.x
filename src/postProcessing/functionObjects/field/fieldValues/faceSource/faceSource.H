/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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
        source          faceZone;   // Type of face source:
                                    // faceZone, patch, sampledSurface
        sourceName      f0;         // faceZone or patch, see below
        operation       sum;
        fields
        (
            p
            phi
            U
        );
    }

    source:
    - faceZone       : requires a 'sourceName' entry to specify the faceZone
    - patch          :                    ""                        patch
    - sampledSurface : requires a 'sampledSurfaceDict' subdictionary. See e.g.
                       sampleDict.

    operation is one of:
      - none
      - sum
      - areaAverage
      - areaIntegrate
      - weightedAverage
      - min
      - max

    Notes:
    - faces on empty patches get ignored
    - if the field is a volField the faceZone can only consist of boundary
    faces.
    - all fields get oriented according to the faceZone (so you might e.g. see
    negative pressure)
    - using sampledSurfaces:
            - they do not do surface fields
            - they use cell values - they do not do any interpolation.
            - take care when using isoSurfaces - these might have duplicate
              triangles so integration might be wrong

SourceFiles
    faceSource.C

\*---------------------------------------------------------------------------*/

#ifndef faceSource_H
#define faceSource_H

#include "NamedEnum.H"
#include "fieldValue.H"
#include "surfaceFieldsFwd.H"
#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class sampledSurface;

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
            stPatch,
            stSampledSurface
        };

        //- Source type names
        static const NamedEnum<sourceType, 3> sourceTypeNames_;


        //- Operation type enumeration
        enum operationType
        {
            opNone,
            opSum,
            opAreaAverage,
            opAreaIntegrate,
            opWeightedAverage,
            opMin,
            opMax
        };

        //- Operation type names
        static const NamedEnum<operationType, 7> operationTypeNames_;


private:

    // Private member functions

        //- Set faces to evaluate based on a face zone
        void setFaceZoneFaces();

        //- Set faces to evaluate based on a patch
        void setPatchFaces();

        //- Set faces according to sampledSurface
        void sampledSurfaceFaces(const dictionary&);


protected:

    // Protected data

        //- Source type
        sourceType source_;

        //- Operation to apply to values
        operationType operation_;

        //- Weight field name - only used for opWeightedAverage mode
        word weightFieldName_;

        //- Global number of faces
        label nFaces_;

        // If operating on mesh faces (faceZone,patch)

            //- Local list of face IDs
            labelList faceId_;

            //- Local list of patch ID per face
            labelList facePatchId_;

            //- List of +1/-1 representing face flip map (1 use as is, -1 negate)
            labelList faceSign_;

        // If operating on sampledSurface

            //- underlying sampledSurface
            autoPtr<sampledSurface> surfacePtr_;


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
            inline const labelList& faceSign() const;


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
                const GeometricField<Type, fvsPatchField, surfaceMesh>& field,
                const bool applyOrientation
            ) const;

            //- Filter a volume field according to faceIds
            template<class Type>
            tmp<Field<Type> > filterField
            (
                const GeometricField<Type, fvPatchField, volMesh>& field,
                const bool applyOrientation
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
