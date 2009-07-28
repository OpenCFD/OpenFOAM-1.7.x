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
    Foam::ensightPart

Description
    Base class for ensightPartCells and ensightPartFaces

SourceFiles
    ensightPart.C
    ensightPartIO.C
    ensightPartI.H

\*---------------------------------------------------------------------------*/

#ifndef ensightPart_H
#define ensightPart_H

#include "ensightFile.H"
#include "ensightGeoFile.H"
#include "typeInfo.H"
#include "labelList.H"
#include "polyMesh.H"
#include "Field.H"
#include "IOPtrList.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ensightPart Declaration
\*---------------------------------------------------------------------------*/

class ensightPart
{
    // Private data

        // Static data members
        static List<word> elemTypes_;


protected:

    // Protected data

        //- part number
        label number_;

        //- part name (or description)
        string name_;

        //- simple labelList with a name
        labelListList elemLists_;

        //- start offset for elemLists_
        label offset_;

        //- number of elements in this part
        label size_;

        //- cell or face data
        bool isCellData_;

        //- material id (numeric)
        label matId_;

        //- mesh reference used
        const polyMesh* meshPtr_;


    // Protected Classes

        //- track the points used by the part and map global to local indices
        class localPoints
        {
        public:
            //- number of points used
            label nPoints;

            //- map global to local indices
            labelList list;

            // null constructor
            localPoints()
            :
                nPoints(0),
                list(0)
            {}

            // construct for mesh points
            localPoints(const polyMesh& pMesh)
            :
                nPoints(0),
                list(pMesh.points().size(), -1)
            {}
        };


    // Protected Member Functions

        //- reconstruct contents from Istream
        void reconstruct(Istream&);

        //- check for fully defined fields
        bool isFieldDefined(const List<scalar>&) const;

        //- write the part header
        void writeHeader(ensightFile&, bool withDescription=false) const;

        //- write a scalar field for idList
        void writeFieldList
        (
            ensightFile& os,
            const List<scalar>& field,
            const List<label>& idList
        ) const;

        //- track points used
        virtual localPoints calcLocalPoints() const
        {
            return localPoints();
        }

        //- write connectivities
        virtual void writeConnectivity
        (
            ensightGeoFile& os,
            const string& key,
            const labelList& idList,
            const labelList& pointMap
        ) const
        {}


public:

    //- Runtime type information
    TypeName("ensightPart");


    // Constructors

        //- Construct null
        ensightPart();

        //- Construct empty part with number and description
        ensightPart(label partNumber, const string& partDescription);

        //- Construct empty part with number and description
        ensightPart
        (
            label partNumber,
            const string& partDescription,
            const polyMesh& pMesh
        );

        //- Construct as copy
        ensightPart(const ensightPart&);


    // Selectors

        // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            ensightPart,
            istream,
            (
                Istream& is
            ),
            (is)
        );

        //- Construct and return clone
        autoPtr<ensightPart> clone() const
        {
            return autoPtr<ensightPart>(new ensightPart(*this));
        };

        //- Construct on freestore from Istream
        static autoPtr<ensightPart> New(Istream& is);


    //- Destructor
    virtual ~ensightPart();


    // Static members

        virtual List<word> const& elementTypes() const
        {
            return elemTypes_;
        }


    // Access

        //- number of elements in this part
        label size() const
        {
            return size_;
        }

        //- represents cell data
        bool isCellData() const
        {
            return isCellData_;
        }

        //- represents face data
        bool isFaceData() const
        {
            return !isCellData_;
        }

        //- part number
        label number() const
        {
            return number_;
        }

        //- part name or description
        const string& name() const
        {
            return name_;
        }

        //- material id
        label materialId() const
        {
            return matId_;
        }

        //- non-const access
        void name(const string& value)
        {
            name_ = value;
        }

        void materialId(const label value)
        {
            matId_ = value;
        }

        //- offset for element ids
        label offset() const
        {
            return offset_;
        }


    // Edit

        //- renumber elements
        void renumber(labelList const&);

        //- write summary information about the object
        bool writeSummary(Ostream&) const;

        //- write reconstruction information for the object
        bool writeData(Ostream&) const;

        //- write geometry
        void writeGeometry(ensightGeoFile&) const;

        //- write scalar field
        void writeScalarField
        (
            ensightFile&,
            const List<scalar>& field
        ) const;

        //- write vector field components
        void writeVectorField
        (
            ensightFile&,
            const List<scalar>& field0,
            const List<scalar>& field1,
            const List<scalar>& field2
        ) const;


        //- write generalized field components
        template <class Type>
        void writeField
        (
            ensightFile&,
            const Field<Type>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ensightPart&)
        {
            notImplemented("ensightPart::operator=(const ensightPart&)");
        }


    // IOstream Operators

        //- write data (reconstruction information)
        friend Ostream& operator<<(Ostream&, const ensightPart&);

        //- write geometry
        friend ensightGeoFile& operator<<
        (
            ensightGeoFile&,
            const ensightPart&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ensightPartI.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
