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
    Foam::autoHexMeshDriver

Description
    main meshing driver.

SourceFiles
    autoHexMeshDriver.C

\*---------------------------------------------------------------------------*/

#ifndef autoHexMeshDriver_H
#define autoHexMeshDriver_H

#include "autoPtr.H"
#include "dictionary.H"
#include "wallPoint.H"
#include "searchableSurfaces.H"
#include "refinementSurfaces.H"
#include "shellSurfaces.H"
#include "meshRefinement.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class fvMesh;

/*---------------------------------------------------------------------------*\
                           Class autoHexMeshDriver Declaration
\*---------------------------------------------------------------------------*/

class autoHexMeshDriver
{
    // Static data members

        //- Extrusion controls
        enum extrudeMode
        {
            NOEXTRUDE,      /*!< Do not extrude. No layers added. */
            EXTRUDE,        /*!< Extrude */
            EXTRUDEREMOVE   /*!< Extrude but afterwards remove added */
                            /*!< faces locally */
        };


    // Private classes

        //- Combine operator class for equalizing displacements.
        class minMagEqOp
        {
        public:

            void operator()(vector& x, const vector& y) const
            {
                if (magSqr(y) < magSqr(x))
                {
                    x = y;
                }
            }
        };

        //- Combine operator class to combine normal with other normal.
        class nomalsCombine
        {
        public:

            void operator()(vector& x, const vector& y) const
            {
                if (y != wallPoint::greatPoint)
                {
                    if (x == wallPoint::greatPoint)
                    {
                        x = y;
                    }
                    else
                    {
                        x *= (x&y);
                    }
                }
            }
        };



    // Private data

        //- Reference to mesh
        fvMesh& mesh_;

        //- Input dictionary
        const dictionary dict_;

        //- Debug level
        const label debug_;

        //- Merge distance
        const scalar mergeDist_;


        //- All surface based geometry
        autoPtr<searchableSurfaces> allGeometryPtr_;

        //- Shells (geometry for inside/outside refinement)
        autoPtr<shellSurfaces> shellsPtr_;

        //- Surfaces (geometry for intersection based refinement)
        autoPtr<refinementSurfaces> surfacesPtr_;

        //- Per refinement surface region the patch
        labelList globalToPatch_;

        //- Mesh refinement engine
        autoPtr<meshRefinement> meshRefinerPtr_;

        //- Decomposition engine
        autoPtr<decompositionMethod> decomposerPtr_;

        //- Mesh distribution engine
        autoPtr<fvMeshDistribute> distributorPtr_;



    // Private Member Functions

        //- Calculate merge distance. Check against writing tolerance.
        scalar getMergeDistance(const scalar mergeTol) const;

        //static void orientOutside(PtrList<searchableSurface>&);

        //- Disallow default bitwise copy construct
        autoHexMeshDriver(const autoHexMeshDriver&);

        //- Disallow default bitwise assignment
        void operator=(const autoHexMeshDriver&);

public:

    //- Runtime type information
    ClassName("autoHexMeshDriver");


    // Constructors

        //- Construct from dictionary and mesh to modify
        autoHexMeshDriver
        (
            fvMesh& mesh,
            const bool overwrite,
            const dictionary& meshDict,
            const dictionary& decomposeDict
        );


    // Member Functions

        // Access

            //- reference to mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }
            fvMesh& mesh()
            {
                return mesh_;
            }

            //- Surfaces to base refinement on
            const refinementSurfaces& surfaces() const
            {
                return surfacesPtr_();
            }

            //- Surfaces to volume refinement on
            const shellSurfaces& shells() const
            {
                return shellsPtr_();
            }

            //- Per refinementsurface, per region the patch
            const labelList& globalToPatch() const
            {
                return globalToPatch_;
            }


        // Meshing

            //- Write mesh
            void writeMesh(const string&) const;

            //- Do all : refine, snap, layers
            void doMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
