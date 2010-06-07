/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::directMappedWallPolyPatch

Description
    Determines a mapping between patch face centres and mesh cell or face
    centres and processors they're on.

Note
    Storage is not optimal. It stores all face centres and cells on all
    processors to keep the addressing calculation simple.

SourceFiles
    directMappedWallPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef directMappedWallPolyPatch_H
#define directMappedWallPolyPatch_H

#include "wallPolyPatch.H"
#include "directMappedPatchBase.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class directMappedWallPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class directMappedWallPolyPatch
:
    public wallPolyPatch,
    public directMappedPatchBase
{

protected:

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry();

        //- Calculate the patch geometry
        virtual void calcGeometry();

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&);

        //- Correct patches after moving points
        virtual void movePoints(const pointField&);

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh();

        //- Update of the patch topology
        virtual void updateMesh();


public:

    //- Runtime type information
    TypeName("directMappedWall");


    // Constructors

        //- Construct from components
        directMappedWallPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from components
        directMappedWallPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const word& sampleRegion,
            const directMappedPatchBase::sampleMode mode,
            const word& samplePatch,
            const vectorField& offset,
            const polyBoundaryMesh& bm
        );

        //- Construct from components. Uniform offset.
        directMappedWallPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const word& sampleRegion,
            const directMappedPatchBase::sampleMode mode,
            const word& samplePatch,
            const vector& offset,
            const polyBoundaryMesh& bm
        );

        //- Construct from dictionary
        directMappedWallPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct as copy, resetting the boundary mesh
        directMappedWallPolyPatch
        (
            const directMappedWallPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        directMappedWallPolyPatch
        (
            const directMappedWallPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>(new directMappedWallPolyPatch(*this, bm));
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new directMappedWallPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~directMappedWallPolyPatch();


    // Member functions

        //- Write the polyPatch data as a dictionary
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
