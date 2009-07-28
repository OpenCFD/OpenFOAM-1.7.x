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
    Foam::fvBoundaryMesh

Description
    Foam::fvBoundaryMesh

SourceFiles
    fvBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef fvBoundaryMesh_H
#define fvBoundaryMesh_H

#include "fvPatchList.H"
#include "lduInterfacePtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                           Class fvBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class fvBoundaryMesh
:
    public fvPatchList
{
private:

    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- Disable default copy construct
        fvBoundaryMesh(const fvBoundaryMesh&);

        //- Disallow assignment
        void operator=(const fvBoundaryMesh&);

        //- Add fvPatches corresponding to the given polyBoundaryMesh
        void addPatches(const polyBoundaryMesh&);


protected:

        //- Update boundary based on new polyBoundaryMesh
        void readUpdate(const polyBoundaryMesh&);


public:

    friend class fvMesh;


    // Constructors

        //- Construct with zero size
        fvBoundaryMesh
        (
            const fvMesh&
        );

        //- Construct from polyBoundaryMesh
        fvBoundaryMesh
        (
            const fvMesh&,
            const polyBoundaryMesh&
        );


    // Member functions

        // Access

            //- Return the mesh reference
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return a list of pointers for each patch
            //  with only those pointing to interfaces being set
            lduInterfacePtrsList interfaces() const;


        //- Correct patches after moving points
        void movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
