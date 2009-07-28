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
    Foam::cuttingPlane

Description
    Constructs plane through mesh.

    No attempt at resolving degenerate cases. Since the cut faces are
    usually quite ugly, they will always be triangulated.

Note
    When the cutting plane coincides with a mesh face, the cell edge on the
    positive side of the plane is taken.

SourceFiles
    cuttingPlane.C

\*---------------------------------------------------------------------------*/

#ifndef cuttingPlane_H
#define cuttingPlane_H

#include "plane.H"
#include "pointField.H"
#include "faceList.H"
#include "MeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class primitiveMesh;

/*---------------------------------------------------------------------------*\
                       Class cuttingPlane Declaration
\*---------------------------------------------------------------------------*/

class cuttingPlane
:
    public plane,
    public MeshedSurface<face>
{
    //- Private typedefs for convenience
        typedef MeshedSurface<face> MeshStorage;

    // Private data

        //- List of cells cut by the plane
        labelList cutCells_;

    // Private Member Functions

        //- Determine cut cells, possibly restricted to a list of cells
        void calcCutCells
        (
            const primitiveMesh&,
            const scalarField& dotProducts,
            const UList<label>& cellIdLabels = UList<label>::null()
        );

        //- Determine intersection points (cutPoints).
        void intersectEdges
        (
            const primitiveMesh&,
            const scalarField& dotProducts,
            List<label>& edgePoint
        );

        //- Walk circumference of cell, starting from startEdgeI crossing
        //  only cut edges. Record cutPoint labels in faceVerts.
        static bool walkCell
        (
            const primitiveMesh&,
            const UList<label>& edgePoint,
            const label cellI,
            const label startEdgeI,
            DynamicList<label>& faceVerts
        );

        //- Determine cuts for all cut cells.
        void walkCellCuts
        (
            const primitiveMesh& mesh,
            const UList<label>& edgePoint
        );


protected:

    // Constructors

        //- Construct plane description without cutting
        cuttingPlane(const plane&);

    // Protected Member Functions

        //- recut mesh with existing planeDesc, restricted to a list of cells
        void reCut
        (
            const primitiveMesh&,
            const UList<label>& cellIdLabels = UList<label>::null()
        );

        //- remap action on triangulation or cleanup
        virtual void remapFaces(const UList<label>& faceMap);

public:

    // Constructors

        //- Construct from plane and mesh reference,
        //  possibly restricted to a list of cells
        cuttingPlane
        (
            const plane&,
            const primitiveMesh&,
            const UList<label>& cellIdLabels = UList<label>::null()
        );


    // Member Functions

        //- Return plane used
        const plane& planeDesc() const
        {
            return static_cast<const plane&>(*this);
        }

        //- Return List of cells cut by the plane
        const labelList& cutCells() const
        {
            return cutCells_;
        }

        //- Return true or false to question: have any cells been cut?
        bool cut() const
        {
            return cutCells_.size();
        }

        //- Sample the cell field
        template<class Type>
        tmp<Field<Type> > sample(const Field<Type>&) const;

        template<class Type>
        tmp<Field<Type> > sample(const tmp<Field<Type> >&) const;


    // Member Operators

        void operator=(const cuttingPlane&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cuttingPlaneTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
