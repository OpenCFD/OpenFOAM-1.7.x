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
    Foam::fvSurfaceMapper

Description
    FV surface mapper.

SourceFiles
    fvSurfaceMapper.C

\*---------------------------------------------------------------------------*/

#ifndef fvSurfaceMapper_H
#define fvSurfaceMapper_H

#include "morphFieldMapper.H"
#include "fvMesh.H"
#include "faceMapper.H"
#include "HashSet.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                         Class fvSurfaceMapper Declaration
\*---------------------------------------------------------------------------*/

class fvSurfaceMapper
:
    public morphFieldMapper
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to face mapper
        const faceMapper& faceMap_;


    // Demand-driven private data

        //- Direct addressing (only one for of addressing is used)
        mutable labelList* directAddrPtr_;

        //- Interpolated addressing (only one for of addressing is used)
        mutable labelListList* interpolationAddrPtr_;

        //- Interpolation weights
        mutable scalarListList* weightsPtr_;

        //- Inserted faces
        mutable labelList* insertedObjectLabelsPtr_;



    // Private Member Functions

        //- Disallow default bitwise copy construct
        fvSurfaceMapper(const fvSurfaceMapper&);

        //- Disallow default bitwise assignment
        void operator=(const fvSurfaceMapper&);


        //- Calculate addressing
        void calcAddressing() const;

        //- Clear out local storage
        void clearOut();


public:

    // Constructors

        //- Construct from components
        fvSurfaceMapper
        (
            const fvMesh& mesh,
            const faceMapper& fMapper
        );


    // Destructor

        virtual ~fvSurfaceMapper();


    // Member Functions

        //- Return size
        virtual label size() const
        {
            return mesh_.nInternalFaces();
        }

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const
        {
            return faceMap_.internalSizeBeforeMapping();
        }

        //- Is the mapping direct
        virtual bool direct() const
        {
            return faceMap_.direct();
        }

        //- Return direct addressing
        virtual const unallocLabelList& directAddressing() const;

        //- Return interpolated addressing
        virtual const labelListList& addressing() const;

        //- Return interpolaion weights
        virtual const scalarListList& weights() const;

        //- Are there any inserted faces
        virtual bool insertedObjects() const
        {
            return faceMap_.insertedObjects();
        }

        //- Return list of inserted faces
        virtual const labelList& insertedObjectLabels() const;

        //- Return flux flip map
        const labelHashSet& flipFaceFlux() const
        {
            return faceMap_.flipFaceFlux();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
