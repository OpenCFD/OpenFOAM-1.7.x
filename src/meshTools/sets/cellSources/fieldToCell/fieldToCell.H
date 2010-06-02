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
    Foam::fieldToCell

Description
    A topoSetSource to select cells based on field values.

SourceFiles
    fieldToCell.C

\*---------------------------------------------------------------------------*/

#ifndef fieldToCell_H
#define fieldToCell_H

#include "topoSetSource.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fieldToCell Declaration
\*---------------------------------------------------------------------------*/

class fieldToCell
:
    public topoSetSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;


        //- Name of volScalarField
        word fieldName_;

        //- min to subset field values with
        scalar min_;

        //- max ,,
        scalar max_;


    // Private Member Functions

        //- Depending on field values add to or delete from cellSet.
        void applyToSet
        (
            const topoSetSource::setAction action,
            const scalarField& field,
            topoSet& set
        ) const;

public:

    //- Runtime type information
    TypeName("fieldToCell");

    // Constructors

        //- Construct from components
        fieldToCell
        (
            const polyMesh& mesh,
            const word& fieldName,
            const scalar min,
            const scalar max
        );

        //- Construct from dictionary
        fieldToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        fieldToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    // Destructor

        virtual ~fieldToCell();


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
