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
    Foam::wallDistReflection

Description
    Wall distance and reflection vector calculation. See wallDist.H

SourceFiles
    wallDistReflection.C

\*---------------------------------------------------------------------------*/

#ifndef wallDistReflection_H
#define wallDistReflection_H

#include "reflectionVectors.H"
#include "wallDistData.H"
#include "wallPointData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class wallDistReflection Declaration
\*---------------------------------------------------------------------------*/

class wallDistReflection
:
    public reflectionVectors,
    public wallDistData<wallPointData<vector> >
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        wallDistReflection(const wallDistReflection&);

        //- Disallow default bitwise assignment
        void operator=(const wallDistReflection&);


public:

    // Constructors

        //- Construct from mesh and flag whether or not to correct wall.
        //  Calculate for all cells. correctWalls : correct wall (face&point)
        //  cells for correct distance, searching neighbours.
        wallDistReflection(const fvMesh& mesh, bool correctWalls = true);


    // Member Functions

        //- Correct for mesh geom/topo changes
        void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
