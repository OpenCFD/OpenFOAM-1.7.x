/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::ensightParts

Description
    A collection of several ensightPart elements

SourceFiles
    ensightParts.C
    ensightPartsI.H

\*---------------------------------------------------------------------------*/

#ifndef ensightParts_H
#define ensightParts_H

#include "ensightPart.H"
#include "ensightPartFaces.H"
#include "ensightPartCells.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ensightParts Declaration
\*---------------------------------------------------------------------------*/

class ensightParts
{
    // Private Data

        //- list of parts
        PtrList<ensightPart> partsList_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        ensightParts(const ensightParts&);

        //- Disallow default bitwise assignment
        void operator=(const ensightParts&);


public:

    // Constructors

        //- Construct from polyMesh
        ensightParts(const polyMesh&);

        //- Construct from IOobject
        ensightParts(const IOobject&);


    //- Destructor
    ~ensightParts();


    // Member functions

        //- clear old information and construct anew from polyMesh
        void recalculate(const polyMesh&);

        //- renumber elements
        void renumber
        (
            const labelList& origCellId,
            const labelList& origFaceId
        );

        //- number of parts
        label size() const
        {
            return partsList_.size();
        }

        //- write the geometry
        void writeGeometry(ensightGeoFile&) const;

        //- write summary information about the objects
        bool writeSummary(Ostream&) const;

        //- write the lists
        void writeData(Ostream&) const;

        //- write scalar field
        void writeScalarField
        (
            ensightFile&,
            const List<scalar>& field,
            bool useFaceData = false
        ) const;

        //- write vector field components
        void writeVectorField
        (
            ensightFile&,
            const List<scalar>& field0,
            const List<scalar>& field1,
            const List<scalar>& field2,
            bool useFaceData = false
        ) const;


        //- write generalized field components
        template <class Type>
        void writeField
        (
            ensightFile&,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;


    // Friend Operators

        friend ensightGeoFile& operator<<
        (
            ensightGeoFile&,
            const ensightParts&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ensightPartsI.H"
#endif

#endif

// ************************************************************************* //
