/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    Foam::TimeActivatedExplicitSource

Description
    Time activated explicit source.

    Sources described by:

    {
        active          true;      // on/off switch
        timeStart       0.2;       // start time
        duration        2.0;       // duration
        selectionMode   points;    // cellSet/cellZone/all
        volumeMode      absolute;  // specific

        fieldData                  // field data - usage for multiple fields
        (
            (H2O 0.005)
        );

        fieldData       0.005;     // field data - usage for single field

        points                     // list of points when selectionMode = points
        (
            (2.75 0.5 0)
        );

        cellSet         c0;        // cellSet name when selectionMode=cellSet
        cellZone        c0;        // cellZone name when selectionMode=cellZone
    }

SourceFiles
    TimeActivatedExplicitSource.C

\*---------------------------------------------------------------------------*/

#ifndef TimeActivatedExplicitSource_H
#define TimeActivatedExplicitSource_H

#include "Tuple2.H"
#include "cellSet.H"
#include "volFieldsFwd.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class fvMesh;

template<class Type>
class TimeActivatedExplicitSource;

// Forward declaration of friend functions

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const TimeActivatedExplicitSource<Type>&
);

/*---------------------------------------------------------------------------*\
               Class TimeActivatedExplicitSource Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TimeActivatedExplicitSource
{
public:

    // Public data

        //- Enumeration for selection mode types
        enum selectionModeType
        {
            smPoints,
            smCellSet,
            smCellZone,
            smAll
        };

        //- Word list of selection mode type names
        static const wordList selectionModeTypeNames_;

        //- Enumeration for volume types
        enum volumeModeType
        {
            vmAbsolute,
            vmSpecific
        };

        //- Word list of volume mode type names
        static const wordList volumeModeTypeNames_;


protected:

    // Protected data

        typedef Tuple2<word, Type> fieldNameValuePair;

        //- Source name
        word name_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Source active flag
        bool active_;

        //- Time start
        scalar timeStart_;

        //- Duration
        scalar duration_;

        //- Volume mode
        volumeModeType volumeMode_;

        //- Cell selection mode
        selectionModeType selectionMode_;

        //- List of points for "points" selectionMode
        List<point> points_;

        //- Name of cell set for "cellSet" and "cellZone" selectionMode
        word cellSetName_;

        //- Set of cells to apply source to
        labelList cells_;

        //- Sum of cell volumes
        scalar V_;

        //- List of source field name vs value pairs
        List<fieldNameValuePair> fieldData_;

        //- Map of fields ids from supplied fields to local field source ids
        labelList fieldIds_;


    // Protected functions

        //- Helper function to convert from a word to a selectionModeType
        selectionModeType wordToSelectionModeType(const word& smtName) const;

        //- Helper function to convert from a word to a volumeModeType
        volumeModeType wordToVolumeModeType(const word& vtName) const;

        //- Helper function to convert from a selectionModeType to a word
        word selectionModeTypeToWord(const selectionModeType& smtType) const;

        //- Helper function to convert from a volumeModeType to a word
        word volumeModeTypeToWord(const volumeModeType& vtType) const;

        //- Set the cellSet or points selection
        void setSelection(const dictionary& dict);

        //- Set the local field data
        void setFieldData(const dictionary& dict, const wordList& fieldNames);

        //- Set the cell set based on the user input selection mode
        void setCellSet();


public:

    // Constructors

        //- Construct from components
        TimeActivatedExplicitSource
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh,
            const wordList& fieldNames
        );

        //- Return clone
        autoPtr<TimeActivatedExplicitSource> clone() const
        {
            notImplemented
            (
                "autoPtr<TimeActivatedExplicitSource> clone() const"
            );
            return autoPtr<TimeActivatedExplicitSource>(NULL);
        }

        //- Return pointer to new TimeActivatedExplicitSource object created
        //  on the freestore from an Istream
        class iNew
        {
            //- Reference to the mesh database
            const fvMesh& mesh_;

            //- List of field names
            const wordList& fieldNames_;


        public:

            iNew
            (
                const fvMesh& mesh,
                const wordList& fieldNames
            )
            :
                mesh_(mesh),
                fieldNames_(fieldNames)
            {}

            autoPtr<TimeActivatedExplicitSource> operator()(Istream& is) const
            {
                const word name(is);
                const dictionary dict(is);

                return autoPtr<TimeActivatedExplicitSource>
                (
                    new TimeActivatedExplicitSource
                    (
                        name,
                        dict,
                        mesh_,
                        fieldNames_
                    )
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to the source name
            inline const word& name() const;

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return const access to the source active flag
            inline bool active() const;

            //- Return const access to the time start
            inline scalar timeStart() const;

            //- Return const access to the duration
            inline scalar duration() const;

            //- Return const access to the time end
            inline scalar timeEnd() const;

            //- Return const access to the volume mode
            inline const volumeModeType& volumeMode() const;

            //- Return const access to the cell selection mode
            inline const selectionModeType& selectionMode() const;

            //- Return const access to the list of points for "points"
            //  selectionMode
            inline const List<point>& points() const;

            //- Return const access to the name of cell set for "cellSet"
            //  selectionMode
            inline const word& cellSetName() const;

            //- Return const access to the total cell volume
            inline scalar V() const;

            //- Return const access to the cell set
            inline const labelList& cells() const;

            //- Return const access to the source field name vs value pairs
            inline const List<fieldNameValuePair>& fieldData() const;

            //- Return const access to the the map of fields ids from supplied
            //  fields to local field source ids
            inline const labelList& fieldIds() const;


        // Edit

            //- Return access to the source name
            inline word& name();

            //- Return access to the source active flag
            inline bool& active();

            //- Return access to the time start
            inline scalar& timeStart();

            //- Return access to the duration
            inline scalar& duration();

            //- Return access to the volume mode
            inline volumeModeType& volumeMode();

            //- Return access to the cell selection mode
            inline selectionModeType& selectionMode();

            //- Return access to the list of points for "points" selectionMode
            inline List<point>& points();

            //- Return access to the name of cell set for "cellSet"
            //  selectionMode
            inline word& cellSetName();

            //- Return access to the total cell volume
            inline scalar& V();

            //- Return access to the cell set
            inline labelList& cells();

            //- Return access to the source field name vs value pairs
            inline List<fieldNameValuePair>& fieldData();

            //- Return access to the the map of fields ids from supplied
            //  fields to local field source ids
            inline labelList& fieldIds();


        // Evaluation

            //- Add the source contribution to field Su
            void addToField
            (
                DimensionedField<Type, volMesh>& Su,
                const label fieldI
            );


        // I-O

            //- Write the source properties
            void writeData(Ostream&) const;

            //- Ostream operator
            friend Ostream& operator<< <Type>
            (
                Ostream& os,
                const TimeActivatedExplicitSource& source
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "TimeActivatedExplicitSource.C"
#   include "TimeActivatedExplicitSourceIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TimeActivatedExplicitSourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
