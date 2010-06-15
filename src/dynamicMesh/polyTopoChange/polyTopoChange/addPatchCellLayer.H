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
    Foam::addPatchCellLayer

Description
    Adds layers of cells to outside of polyPatch.

    Call setRefinement with offset vector for every patch point and number
    of layers per patch face and number of layers per patch point.
    - offset vector should be zero for any non-manifold point and synchronised
      on coupled points before calling this.
    - offset vector of zero will not add any points.
    - gets supplied the number of extruded layers both per face and per
      point. Usually the point nlayers is the max of surrounding face nlayers.

      point nlayers:
       -  0 : no extrusion. Any surrounding face being extruded becomes 'prism'
       - >0 : should be max of surrounding face nlayers.

    - differing face nlayers: 'termination' : (e.g. from 2 to 4 layers) match
      at original patch face side.

        E.g. 2 boundary faces on patches a,b. 2 layers for a, 3 for b.

    @verbatim
        Was:

           a      b         <- patch of boundary face
        +------+------+
        |      |      |     <- original cells
        +------+------+

        Becomes:

           a      b         <- patch of boundary face
        +------+------+
        +      +------+
        +------+------+
        +------+------+
        |      |      |     <- original cells
        +------+------+
    @endverbatim


    - added faces get same patchID as face they are extruded from
    - 'side' faces (i.e. on the edge of pp) get the patchID of the
    other patch they are connected to.


    E.g. 3 boundary faces on patches a,b. b gets extruded, a doesn't.

    @verbatim
           a      b      b          <- patch of boundary face
        +------+------+------+
        |      |      |      |      <- cells
        +------+------+------+


               ^      ^             <- wanted extrusion vector (none at far right)
           a   |  b   |  b          <- patch of boundary face
        +------+------+------+
        |      |      |      |      <- cells
        +------+------+------+

                  b
               +------+\ b          1. prism cell added onto second b face since
           a  a|      | ----\          only one side gets extruded.
        +------+------+------+      2. side-face gets patch a, not b.
        |      |      |      |
        +------+------+------+
    @endverbatim


SourceFiles
    addPatchCellLayer.C

\*---------------------------------------------------------------------------*/

#ifndef addPatchCellLayer_H
#define addPatchCellLayer_H

#include "labelList.H"
#include "typeInfo.H"
#include "labelPair.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class mapPolyMesh;
class primitiveMesh;
class globalIndex;

/*---------------------------------------------------------------------------*\
                           Class addPatchCellLayer Declaration
\*---------------------------------------------------------------------------*/

class addPatchCellLayer
{
    // Private classes

        // To combineReduce a labelList. Filters out duplicates.
        class uniqueEqOp
        {

        public:

            void operator()(labelList& x, const labelList& y) const
            {
                if (x.empty())
                {
                    if (y.size())
                    {
                        x = y;
                    }
                }
                else
                {
                    forAll(y, yi)
                    {
                        if (findIndex(x, y[yi]) == -1)
                        {
                            label sz = x.size();
                            x.setSize(sz+1);
                            x[sz] = y[yi];
                        }
                    }
                }
            }
        };



    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- For all patchpoints: list of added points (size 0 or nLayers)
        //  First point in list is one nearest to original point in patch,
        //  last one is the new point on the surface.
        labelListList addedPoints_;

        //- For all patchfaces: list of layer faces.
        //  - empty if no face extruded
        //  - first face is original boundary face
        //  - last one is new boundary face.
        labelListList layerFaces_;


    // Private Member Functions

        //- Per patch edge the pp faces (in global indices) using it. Uses
        //  uniqueEqOp() to remove duplicates.
        labelListList calcGlobalEdgeFaces
        (
            const polyMesh& mesh,
            const globalIndex& globalFaces,
            const indirectPrimitivePatch& pp,
            const labelList& meshEdges
        );

        //- Get the face on the other side of the edge.
        static label nbrFace
        (
            const labelListList& edgeFaces,
            const label edgeI,
            const label faceI
        );

        //- Add vertex to face if unique.
        static void addVertex(const label, face&, label& fp);

        bool sameEdgeNeighbour
        (
            const indirectPrimitivePatch& pp,
            const labelListList& globalEdgeFaces,
            const boolList& doneEdge,
            const label thisGlobalFaceI,
            const label nbrGlobalFaceI,
            const label edgeI
        ) const;

        labelPair getEdgeString
        (
            const indirectPrimitivePatch& pp,
            const labelListList& globalEdgeFaces,
            const boolList& doneEdge,
            const label patchFaceI,
            const label globalFaceI
        ) const;


        //- Add face between layer-1 and layer.
        label addSideFace
        (
            const indirectPrimitivePatch&,
            const labelList& patchID,
            const labelListList& addedCells,
            const face& newFace,
            const label ownFaceI,
            const label nbrFaceI,
            const label patchEdgeI,
            const label meshEdgeI,
            const label layerI,
            const label numEdgeFaces,
            const labelList& meshFaces,
            polyTopoChange&
        ) const;


        //- Disallow default bitwise copy construct
        addPatchCellLayer(const addPatchCellLayer&);

        //- Disallow default bitwise assignment
        void operator=(const addPatchCellLayer&);


public:

    //- Runtime type information
    ClassName("addPatchCellLayer");


    // Constructors

        //- Construct from mesh.
        addPatchCellLayer(const polyMesh& mesh);


    // Member Functions


        // Access

            //- Added points per patch point.
            const labelListList& addedPoints() const
            {
                return addedPoints_;
            }

            //- Layer faces per patch face. See above.
            const labelListList& layerFaces() const
            {
                return layerFaces_;
            }

            //- Helper: get added cells per patch face.
            //  addedCells[patchFace] is list of cells added. Last element is
            //  the top cells (i.e. the boundary cell)
            static labelListList addedCells
            (
                const polyMesh&,
                const labelListList& layerFaces
            );

            //- added cells given current mesh & layerfaces.
            labelListList addedCells() const;


        // Edit

            //- Play commands into polyTopoChange to create layers on top
            //  of indirectPrimitivePatch (have to be outside faces).
            //  Gets displacement per patch point.
            //  - nPointLayers : number of layers per (patch)point
            //  - nFaceLayers : number of layers per (patch) face
            //  - firstDisplacement : displacement per point for first
            //    layer of points (i.e. nearest to original mesh). If zero
            //    do not add point.
            //  Layer thicknesses are calculated to constant geometric
            //  expansion. Use expansionRatio 1 for constant size.
            //  Sets addedPoints_ which is per pp point a list of points
            //  added.
            //  Note: firstDisplacement has to be parallel synchronised before
            //        calling this routine. Only if all procs sharing a point
            //        get a cell should firstDisplacement be <> 0
            //  Note: cells get added from owner cells of patch faces
            //        (instead of e.g. from patch faces)
           void setRefinement
           (
               const scalarField& expansionRatio,
               const indirectPrimitivePatch& pp,
               const labelList& nFaceLayers,
               const labelList& nPointLayers,
               const vectorField& firstLayerDisp,
               polyTopoChange& meshMod
            );


            //- Add with constant expansion ratio and same nLayers everywhere
            void setRefinement
            (
                const label nLayers,
                const indirectPrimitivePatch& pp,
                const vectorField& overallDisplacement,
                polyTopoChange& meshMod
            )
            {
                setRefinement
                (
                    scalarField(pp.nPoints(), 1.0), // expansion ration
                    pp,
                    labelList(pp.size(), nLayers),
                    labelList(pp.nPoints(), nLayers),
                    overallDisplacement / nLayers,
                    meshMod
                );
            }


            //- Update any locally stored mesh information. Gets additional
            //  map from new to old patch (since patch needs to be
            //  recreated since has to be on outside).
            void updateMesh
            (
                const mapPolyMesh&,
                const labelList& faceMap,   // new to old patch faces
                const labelList& pointMap   // new to old patch points
            );

        // Helper

            //- Per patch edge the corresponding mesh edge
            static labelList calcMeshEdges
            (
                const primitiveMesh& mesh,
                const indirectPrimitivePatch&
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
