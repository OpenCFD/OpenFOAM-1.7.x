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
    Foam::CFCCellToCellStencil

Description

SourceFiles
    CFCCellToCellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef CFCCellToCellStencil_H
#define CFCCellToCellStencil_H

#include "cellToCellStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class CFCCellToCellStencil Declaration
\*---------------------------------------------------------------------------*/

class CFCCellToCellStencil
:
    public cellToCellStencil
{
    // Private Member Functions

        void calcFaceBoundaryData(labelList& neiGlobal) const;

        void calcCellStencil(labelListList& globalCellCells) const;

        //- Disallow default bitwise copy construct
        CFCCellToCellStencil(const CFCCellToCellStencil&);

        //- Disallow default bitwise assignment
        void operator=(const CFCCellToCellStencil&);

public:

    // Constructors

        //- Construct from mesh
        explicit CFCCellToCellStencil(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
