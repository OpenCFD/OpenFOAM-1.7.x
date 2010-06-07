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
    Foam::globalPointPatch

Description
    Foam::globalPointPatch

SourceFiles
    globalPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef globalPointPatch_H
#define globalPointPatch_H

#include "pointPatch.H"
#include "coupledPointPatch.H"
#include "globalMeshData.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class globalPointPatch Declaration
\*---------------------------------------------------------------------------*/

class globalPointPatch
:
    public pointPatch,
    public coupledPointPatch
{
    // Private data

        // //- Reference to the basic globalMeshData
        // const globalMeshData& globalMeshData_;

        //- Index in the boundary mesh
        label index_;


    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry()
        {}

        //- Calculate the patch geometry
        virtual void calcGeometry()
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&)
        {}

        //- Correct patches after moving points
        virtual void movePoints(const pointField&)
        {}

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh()
        {}

        //- Update of the patch topology
        virtual void updateMesh()
        {}


    // Private Member Functions

        //- Disallow default construct as copy
        globalPointPatch
        (
            const globalPointPatch&
        );

        //- Disallow default assignment
        void operator=(const globalPointPatch&);


public:

    //- Runtime type information
    TypeName("global");

    // Constructors

        //- Construct from components
        globalPointPatch
        (
            const pointBoundaryMesh&,
            const label index
        );


    // Destructor

        virtual ~globalPointPatch();


    // Member functions

        //- Return name
        virtual const word& name() const
        {
            // There can only be a single patch of this type - therefore
            // its name is hard-coded.
            return type();
        }

        //- Return size
        virtual label size() const
        {
            return meshPoints().size();
        }

        //- Return true if running parallel
        virtual bool coupled() const
        {
            if (Pstream::parRun())
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        //- Return number of faces
        virtual label nFaces() const
        {
            return 0;
        }

        //- Return the index of this patch in the pointBoundaryMesh
        virtual label index() const
        {
            return index_;
        }

        //- Return mesh points
        virtual const labelList& meshPoints() const
        {
            return boundaryMesh().mesh().globalData().sharedPointLabels();
        }

        //- Return local points.  Not implemented
        virtual const pointField& localPoints() const
        {
            notImplemented("globalPointPatch::localPoints() const");
            return pointField::null();
        }

        //- Return point normals.  Not implemented
        virtual const vectorField& pointNormals() const
        {
            notImplemented("globalPointPatch::pointNormals() const");
            return vectorField::null();
        }

        //- Return total number of shared points
        virtual label globalPointSize() const
        {
            return boundaryMesh().mesh().globalData().nGlobalPoints();
        }

        //- Return addressing into the global points list
        const labelList& sharedPointAddr() const
        {
            return boundaryMesh().mesh().globalData().sharedPointAddr();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
