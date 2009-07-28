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
    Foam::dynamicFvMesh

Description
    Abstract base class for geometry and/or topology changing fvMesh.

SourceFiles
    dynamicFvMesh.C
    newDynamicFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicFvMesh_H
#define dynamicFvMesh_H

#include "fvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dynamicFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicFvMesh
:
    public fvMesh
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicFvMesh(const dynamicFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicFvMesh");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            dynamicFvMesh,
            IOobject,
            (const IOobject& io),
            (io)
        );


    // Constructors

        //- Construct from objectRegistry, and read/write options
        explicit dynamicFvMesh(const IOobject& io);


    // Selectors

        //- Select null constructed
        static autoPtr<dynamicFvMesh> New(const IOobject& io);


    // Destructor

        virtual ~dynamicFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
