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
    Foam::pointBoundaryMesh

Description
    Foam::pointBoundaryMesh

SourceFiles
    pointBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef pointBoundaryMesh_H
#define pointBoundaryMesh_H

#include "pointPatchList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class pointMesh;
class polyBoundaryMesh;
class globalPointPatch;

/*---------------------------------------------------------------------------*\
                       Class pointBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class pointBoundaryMesh
:
    public pointPatchList
{
    // Private data

        //- Reference to mesh
        const pointMesh& mesh_;


    // Private Member Functions

        //- Calculate the geometry for the patches (transformation tensors etc.)
        void calcGeometry();

        //- Disallow default bitwise copy construct
        pointBoundaryMesh(const pointBoundaryMesh&);

        //- Disallow default bitwise assignment
        void operator=(const pointBoundaryMesh&);


public:

    //- Declare friendship with pointMesh
    friend class pointMesh;


    // Constructors

        //- Construct from polyBoundaryMesh
        pointBoundaryMesh
        (
            const pointMesh&,
            const polyBoundaryMesh&
        );


    // Member functions

        //- Return the mesh reference
        const pointMesh& mesh() const
        {
            return mesh_;
        }

        //- Return reference to globalPointPatch
        const globalPointPatch& globalPatch() const;

        //- Correct polyBoundaryMesh after moving points
        void movePoints(const pointField&);

        //- Correct polyBoundaryMesh after topology update
        void updateMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
