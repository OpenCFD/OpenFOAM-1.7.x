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
    Foam::centredCFCFaceToCellStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredCFCFaceToCellStencilObject_H
#define centredCFCFaceToCellStencilObject_H

#include "extendedCentredFaceToCellStencil.H"
#include "CFCFaceToCellStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class centredCFCFaceToCellStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredCFCFaceToCellStencilObject
:
    public MeshObject<fvMesh, centredCFCFaceToCellStencilObject>,
    public extendedCentredFaceToCellStencil
{

public:

    TypeName("centredCFCFaceToCellStencil");

    // Constructors

        //- Construct from uncompacted face stencil
        explicit centredCFCFaceToCellStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject<fvMesh, centredCFCFaceToCellStencilObject>(mesh),
            extendedCentredFaceToCellStencil(CFCFaceToCellStencil(mesh))
        {}


    //- Destructor
    virtual ~centredCFCFaceToCellStencilObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
