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
    Foam::coupledPolyPatch

Description
    The coupledPolyPatch is an abstract base class for patches that couple
    regions of the computational domain e.g. cyclic and processor-processor
    links.

SourceFiles
    coupledPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef coupledPolyPatch_H
#define coupledPolyPatch_H

#include "polyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class coupledPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledPolyPatch
:
    public polyPatch
{
    // Private data

        //- offset (distance) vector from one side of the couple to the other
        mutable vectorField separation_;

        //- Face transformation tensor
        mutable tensorField forwardT_;

        //- Neighbour-cell transformation tensor
        mutable tensorField reverseT_;

public:

    // Static data members

        //- Relative tolerance (for geometric matching).
        static scalar matchTol;


protected:

    // Protected Member Functions

        //- Calculate the transformation tensors
        //  smallDist : matching distance per face
        //  absTol    : absolute error in normal
        void calcTransformTensors
        (
            const vectorField& Cf,
            const vectorField& Cr,
            const vectorField& nf,
            const vectorField& nr,
            const scalarField& smallDist,
            const scalar absTol = matchTol
        ) const;

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry() = 0;

        //- Calculate the patch geometry
        virtual void calcGeometry() = 0;

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&) = 0;

        //- Correct patches after moving points
        virtual void movePoints(const pointField&) = 0;

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh() = 0;

        //- Update of the patch topology
        virtual void updateMesh() = 0;


        //- Write point in OBJ format
        static void writeOBJ(Ostream& os, const point& pt);

        //- Write selected points in OBJ format
        static void writeOBJ(Ostream&, const pointField&, const labelList&);

        //- Write patch
        static void writeOBJ
        (
            const fileName&,
            const UList<face>&,
            const pointField&
        );

        //- Write edge in OBJ format
        static void writeOBJ
        (
            Ostream& os,
            const point& p0,
            const point& p1,
            label& vertI
        );

        //- Calculate face centres
        static pointField calcFaceCentres
        (
            const UList<face>&,
            const pointField&
        );

        //- Get f[0] for all faces
        static pointField getAnchorPoints
        (
            const UList<face>&,
            const pointField&
        );

        //- Is face (in old face labels) in current patch?
        bool inPatch
        (
            const labelList& oldToNew,
            const label oldFaceI
        ) const;

        //- Given list of starts of patches and a face label determine
        //  the patch.
        static label whichPatch
        (
            const labelList& patchStarts,
            const label faceI
        );

        //- Calculate typical tolerance per face. Is currently max distance
        //  from face centre to any of the face vertices.
        static scalarField calcFaceTol
        (
            const UList<face>& faces,
            const pointField& points,
            const pointField& faceCentres
        );

        //- Get the number of vertices face f needs to be rotated such that
        //  its f[0] point aligns with given anchor (within tol).
        static label getRotation
        (
            const pointField& points,
            const face& f,
            const point& anchor,
            const scalar tol
        );


public:

    //- Runtime type information
    TypeName("coupled");


    // Constructors

        //- Construct from components
        coupledPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from dictionary
        coupledPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct as copy, resetting the boundary mesh
        coupledPolyPatch(const coupledPolyPatch&, const polyBoundaryMesh&);

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        coupledPolyPatch
        (
            const coupledPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );


    // Destructor

        virtual ~coupledPolyPatch();


    // Member Functions

        // Access

            //- Return true because this patch is coupled
            virtual bool coupled() const
            {
                return true;
            }


            //- Are the coupled planes separated
            bool separated() const
            {
                return separation_.size();
            }

            //- Return the offset (distance) vector from one side of the couple
            //  to the other
            const vectorField& separation() const
            {
                return separation_;
            }


            //- Are the cyclic planes parallel
            bool parallel() const
            {
                return forwardT_.empty();
            }

            //- Return face transformation tensor
            const tensorField& forwardT() const
            {
                return forwardT_;
            }

            //- Return neighbour-cell transformation tensor
            const tensorField& reverseT() const
            {
                return reverseT_;
            }


        //- Initialize ordering for primitivePatch. Does not
        //  refer to *this (except for name() and type() etc.)
        virtual void initOrder(const primitivePatch&) const = 0;

        //- Return new ordering for primitivePatch.
        //  Ordering is -faceMap: for every face
        //  index of the new face -rotation:for every new face the clockwise
        //  shift of the original face. Return false if nothing changes
        //  (faceMap is identity, rotation is 0), true otherwise.
        virtual bool order
        (
            const primitivePatch&,
            labelList& faceMap,
            labelList& rotation
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
