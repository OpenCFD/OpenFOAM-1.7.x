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
    Foam::thresholdCellFaces

Description
    Selects the mesh cell faces specified by a threshold value.
    Non-triangulated by default.

SourceFiles
    thresholdCellFaces.C

\*---------------------------------------------------------------------------*/

#ifndef thresholdCellFaces_H
#define thresholdCellFaces_H

#include "MeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class thresholdCellFaces Declaration
\*---------------------------------------------------------------------------*/

class thresholdCellFaces
:
    public MeshedSurface<face>
{
    //- Private typedefs for convenience
        typedef MeshedSurface<face> MeshStorage;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- For every face the original cell in mesh
        labelList meshCells_;

    // Private Member Functions

        void calculate
        (
            const scalarField&, 
            const scalar lowerThreshold,
            const scalar upperThreshold,
            const bool triangulate
        );

public:

    //- Runtime type information
    TypeName("thresholdCellFaces");


    // Constructors

        //- Construct from mesh, field and threshold value
        thresholdCellFaces
        (
            const polyMesh&,
            const scalarField&,
            const scalar lowerThreshold,
            const scalar upperThreshold,
            const bool triangulate = false
        );


    // Member Functions

        //- For every face original cell in mesh
        labelList& meshCells()
        {
            return meshCells_;
        }

        //- For every face original cell in mesh
        const labelList& meshCells() const
        {
            return meshCells_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
