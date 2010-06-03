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
    Foam::regionToCell

Description
    TopoSetSource. Select cells belonging to topological connected region
    (that contains given point)

SourceFiles
    regionToCell.C

\*---------------------------------------------------------------------------*/

#ifndef regionToCell_H
#define regionToCell_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class regionToCell Declaration
\*---------------------------------------------------------------------------*/

class regionToCell
:
    public topoSetSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of cellSet to keep to
        word setName_;

        //- Coordinate that is inside connected region
        point insidePoint_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("regionToCell");

    // Constructors

        //- Construct from components
        regionToCell
        (
            const polyMesh& mesh,
            const word& setName,
            const point& insidePoint
        );

        //- Construct from dictionary
        regionToCell
        (
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from Istream
        regionToCell
        (
            const polyMesh& mesh,
            Istream&
        );


    // Destructor

        virtual ~regionToCell();


    // Member Functions

        virtual void applyToSet(const topoSetSource::setAction action, topoSet&)
         const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
