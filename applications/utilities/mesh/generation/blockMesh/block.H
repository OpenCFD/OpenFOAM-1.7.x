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
    Foam::block

Description
    creates a single block of cells from point coordinates,
    numbers of cells in each direction and expansion ratio

SourceFiles
    block.C

\*---------------------------------------------------------------------------*/

#ifndef block_H
#define block_H

#include "pointField.H"
#include "labelList.H"

#include "blockDescriptor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                           Class block Declaration
\*---------------------------------------------------------------------------*/

class block
{
    // Private data

        //- block definition
        blockDescriptor blockDef_;

        //- list of vertices
        pointField vertices_;

        //- list of cells
        labelListList cells_;

        //- boundary patches
        labelListListList boundaryPatches_;


    // private member functions

        label vtxLabel(label i, label j, label k);

        void blockPoints();

        void blockCells();

        void blockBoundary();


public:

    // Constructors

        //- Construct from the block definition
        block(const blockDescriptor&);

        //- Clone
        autoPtr<block> clone() const
        {
            notImplemented("block::clone()");
            return autoPtr<block>(NULL);
        }


    // Member Functions

        // Access

            const blockDescriptor& blockDef() const;
            const pointField& points() const;
            const labelListList& cells() const;
            const labelListListList& boundaryPatches() const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const block&);
};


inline Istream& operator>>(Istream& is, block*)
{
    notImplemented("Istream& operator>>(Istream& is, block*)");
    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
