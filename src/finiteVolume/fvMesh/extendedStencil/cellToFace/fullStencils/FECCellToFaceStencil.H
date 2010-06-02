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
    Foam::FECCellToFaceStencil

Description
    All cells connected via edge to face.

SourceFiles
    FECCellToFaceStencil.C

\*---------------------------------------------------------------------------*/

#ifndef FECCellToFaceStencil_H
#define FECCellToFaceStencil_H

#include "cellToFaceStencil.H"
#include "EdgeMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FECCellToFaceStencil Declaration
\*---------------------------------------------------------------------------*/

class FECCellToFaceStencil
:
    public cellToFaceStencil
{
    // Private Member Functions

        //- Calculates per edge the neighbour data (= edgeCells)
        void calcEdgeBoundaryData
        (
            const boolList& isValidBFace,
            const labelList& boundaryEdges,
            EdgeMap<labelList>& neiGlobal
        ) const;

        void calcFaceStencil(labelListList& faceStencil) const;

public:

    // Constructors

        //- Construct from all cells and boundary faces
        explicit FECCellToFaceStencil(const polyMesh&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
