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
    Foam::basicSource

Description
    Basic source abtract class

    Sources described by:

    source1
    {
        typeModel       actuationDiskSource; // explicitSource
        active          on;            // on/off switch
        timeStart       0.0;           // start time
        duration        1000.0;        // duration
        selectionMode   cellSet;       // cellSet // points //cellZone
        cellSet         c0;            // cellSet name

        actuationDiskSourceCoeffs
        {
            diskDir     (-1 0 0); // orientation of the disk
            Cp          0.53;     // Cp
            Ct          0.58;     // Ct
            diskArea    40;       // disk area
        }
    }

    source2
    {
        typeModel       explicitSource;
        active          on;
        timeStart       0.0;
        duration        1000.0;
        selectionMode   points;
        cellSet         c0;

        explicitSourceCoeffs
        {
            points            // list of points when selectionMode = points
            (
                (-0.088 0.007 -0.02)
                (-0.028 0.007 -0.02)
            );
            volumeMode      specific;  //absolute
            fieldData                  //field data
            {
                k   30.7;
                epsilon  1.5;
            }
        }
    }

SourceFiles
    basicSource.C
    basicSourceIO.C

\*---------------------------------------------------------------------------*/

#ifndef basicSource_H
#define basicSource_H

#include "fvMatrices.H"
#include "cellSet.H"
#include "volFieldsFwd.H"
#include "DimensionedField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class basicSource Declaration
\*---------------------------------------------------------------------------*/

class basicSource
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


protected:

    // Protected data

        //- Source name
        word name_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Dictionary containing the data of the source
        const dictionary& dict_;

        //- Source active flag
        bool active_;

        //- Time start
        scalar timeStart_;

        //- Duration
        scalar duration_;

        //- Cell selection mode
        selectionModeType selectionMode_;

        //- Name of cell set for "cellSet" and "cellZone" selectionMode
        word cellSetName_;

        //- Set of cells to apply source to
        labelList cells_;

        //- Sum of cell volumes
        scalar V_;


    // Protected functions

        //- Helper function to convert from a word to a selectionModeType
        selectionModeType wordToSelectionModeType(const word& smtName) const;

        //- Helper function to convert from a selectionModeType to a word
        word selectionModeTypeToWord(const selectionModeType& smtType) const;

        //- Set the cellSet or points selection
        void setSelection(const dictionary& dict);

        //- Set the cell set based on the user input selection mode
        void setCellSet();


public:

    //- Runtime type information
    TypeName("basicSource");


     // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            basicSource,
            dictionary,
            (
                const word& name,
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (name, dict, mesh)
        );


    // Constructors

        //- Construct from components
        basicSource
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Return clone
        autoPtr<basicSource> clone() const
        {
            notImplemented
            (
                "autoPtr<basicSource> clone() const"
            );
            return autoPtr<basicSource>(NULL);
        }

        //- Return pointer to new basicSource object created
        //  on the freestore from an Istream
        class iNew
        {
            //- Reference to the mesh database
            const fvMesh& mesh_;
            const word& name_;

        public:

            iNew
            (
                const fvMesh& mesh,
                const word& name
            )
            :
                mesh_(mesh),
                name_(name)
            {}

            autoPtr<basicSource> operator()(Istream& is) const
            {
                //const word name(is);
                const dictionary dict(is);

                return autoPtr<basicSource>
                (
                    basicSource::New
                    (
                        name_,
                        dict,
                        mesh_
                    )
                );
            }
        };


    // Selectors

        //- Return a reference to the selected basicSource model
        static autoPtr<basicSource> New
        (
            const word& name,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~basicSource()
    {}


    // Member Functions

        // Access

            //- Return const access to the source name
            inline const word& name() const;

            //- Return const access to the mesh database
            inline const fvMesh& mesh() const;

            //- Return dictionay
            inline const dictionary& dictCoeffs() const;

            //- Return const access to the source active flag
            inline bool active() const;

            //- Return const access to the time start
            inline scalar timeStart() const;

            //- Return const access to the duration
            inline scalar duration() const;

            //- Return const access to the time end
            inline scalar timeEnd() const;

            //- Return const access to the cell selection mode
            inline const selectionModeType& selectionMode() const;

            //- Return const access to the name of cell set for "cellSet"
            //  selectionMode
            inline const word& cellSetName() const;

            //- Return const access to the total cell volume
            inline scalar V() const;

            //- Return const access to the cell set
            inline const labelList& cells() const;


        // Edit

            //- Return access to the source name
            inline word& name();

            //- Return access to the source active flag
            inline bool& active();

            //- Return access to the time start
            inline scalar& timeStart();

            //- Return access to the duration
            inline scalar& duration();

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


        // Checks

            //- Is the source active?
            bool isActive();


        // Evaluation

            //- Add all explicit sources
            virtual void addExplicitSources() = 0;

            //- Add source to scalar field
            virtual void addSu(DimensionedField<scalar, volMesh>& field) = 0;

            //- Add source to vector field
            virtual void addSu(DimensionedField<vector, volMesh>& field) = 0;

            //- Add source term to vector fvMatrix
            virtual void addSu(fvMatrix<vector>& Eqn) = 0;

            //- Add source term to scalar fvMatrix
            virtual void addSu(fvMatrix<scalar>& Eqn) = 0;


        // I-O

            //- Write the source properties
            virtual void writeData(Ostream&) const = 0;

            //- Read source dictionary
            virtual bool read(const dictionary& dict) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicSourceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
