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
    Foam::wallDist

Description
    Calculation of distance to nearest wall for all cells and boundary.
    Uses meshWave to do actual calculation.

    Distance correction:

    if correctWalls = true:
    For each cell with face on wall calculate the true nearest point
    (by triangle decomposition) on that face and do that same for that face's
    pointNeighbours. This will find the true nearest distance in almost all
    cases. Only very skewed cells or cells close to another wall might be
    missed.

    For each cell with only point on wall the same is done except now it takes
    the pointFaces() of the wall point to look for the nearest point.

Note

    correct() : for now does complete recalculation. (which usually is
    ok since mesh is smoothed). However for topology change where geometry
    in most of domain does not change you could think of starting from the
    old cell values. Tried but not done since:
    - meshWave would have to be called with old cellInfo.
      This is List\<wallInfo\> of nCells.
    - cannot construct from distance (y_) only since we don't know a value
      for origin_. (origin_ = GREAT already used to denote illegal value.)
    - so we would have to store a List\<wallInfo\> which unfortunately does
      not get resized/mapped automatically upon mesh changes.

SourceFiles
    wallDist.C

\*---------------------------------------------------------------------------*/

#ifndef wallDist_H
#define wallDist_H

#include "volFields.H"
#include "cellDistFuncs.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class wallDist Declaration
\*---------------------------------------------------------------------------*/

class wallDist
:
    public volScalarField,
    public cellDistFuncs
{


private:

    // Private Member Data

        //- Do accurate distance calculation for near-wall cells.
        bool correctWalls_;

        //- Number of unset cells and faces.
        label nUnset_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        wallDist(const wallDist&);

        //- Disallow default bitwise assignment
        void operator=(const wallDist&);


public:

    // Constructors

        //- Construct from mesh and flag whether or not to correct wall.
        //  Calculate for all cells. correctWalls : correct wall (face&point)
        //  cells for correct distance, searching neighbours.
        wallDist(const fvMesh& mesh, bool correctWalls = true);


    // Destructor

        virtual ~wallDist();


    // Member Functions

        const volScalarField& y() const
        {
            return *this;
        }

        label nUnset() const
        {
            return nUnset_;
        }

        //- Correct for mesh geom/topo changes
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
