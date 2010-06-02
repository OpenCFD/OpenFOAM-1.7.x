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
    Foam::nearWallDistNoSearch

Description
    Distance calculation for cells with face on a wall. Does not search
    anything, just takes normal component of distance.

SourceFiles
    nearWallDistNoSearch.C

\*---------------------------------------------------------------------------*/

#ifndef nearWallDistNoSearch_H
#define nearWallDistNoSearch_H

#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class nearWallDistNoSearch Declaration
\*---------------------------------------------------------------------------*/

class nearWallDistNoSearch
:
    public volScalarField::GeometricBoundaryField
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- Do all calculations.
        void doAll();

        //- Disallow default bitwise copy construct
        nearWallDistNoSearch(const nearWallDistNoSearch&);

        //- Disallow default bitwise assignment
        void operator=(const nearWallDistNoSearch&);


public:

    // Constructors

        //- Construct from components
        nearWallDistNoSearch(const fvMesh& mesh);


    // Destructor

        virtual ~nearWallDistNoSearch();


    // Member Functions

        //- Correct for mesh geom/topo changes
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
