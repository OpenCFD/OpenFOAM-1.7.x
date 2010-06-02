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
    Foam::autoLayerDriver

Description
    All to do with adding layers

SourceFiles
    autoLayerDriver.C

\*---------------------------------------------------------------------------*/

#ifndef autoLayerDriver_H
#define autoLayerDriver_H

#include "meshRefinement.H"
#include "wallPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class removePoints;
class pointSet;
class motionSmoother;
class addPatchCellLayer;
class pointData;
class wallPoint;
class faceSet;
class layerParameters;

/*---------------------------------------------------------------------------*\
                           Class autoLayerDriver Declaration
\*---------------------------------------------------------------------------*/

class autoLayerDriver
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

        //- Mesh+surface
        meshRefinement& meshRefiner_;



    // Private Member Functions


        // Face merging

            //- Merge patch faces. Undo until no checkMesh errors.
            label mergePatchFacesUndo
            (
                const scalar minCos,
                const scalar concaveCos,
                const dictionary&
            );

            //- Remove points.
            autoPtr<mapPolyMesh> doRemovePoints
            (
                removePoints& pointRemover,
                const boolList& pointCanBeDeleted
            );

            //- Restore faces (which contain removed points)
            autoPtr<mapPolyMesh> doRestorePoints
            (
                removePoints& pointRemover,
                const labelList& facesToRestore
            );

            //- Return candidateFaces that are also in set.
            labelList collectFaces
            (
                const labelList& candidateFaces,
                const labelHashSet& set
            ) const;

            //- Pick up faces of cells of faces in set.
            labelList growFaceCellFace(const labelHashSet&) const;

            //- Remove points not used by any face or points used by only
            //  two faces where the edges are in line
            label mergeEdgesUndo(const scalar minCos, const dictionary&);


        // Layers

            //- For debugging: Dump displacement to .obj files
            static void dumpDisplacement
            (
                const fileName&,
                const indirectPrimitivePatch&,
                const vectorField&,
                const List<extrudeMode>&
            );

            //- Check that primitivePatch is not multiply connected.
            //  Collect non-manifold points in pointSet.
            static void checkManifold
            (
                const indirectPrimitivePatch&,
                pointSet& nonManifoldPoints
            );

            //- Check that mesh outside is not multiply connected.
            void checkMeshManifold() const;


            // Static extrusion setup

                //- Unset extrusion on point. Returns true if anything unset.
                static bool unmarkExtrusion
                (
                    const label patchPointI,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Unset extrusion on face. Returns true if anything unset.
                static bool unmarkExtrusion
                (
                    const face& localFace,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- No extrusion at non-manifold points.
                void handleNonManifolds
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion on feature edges. Assumes non-manifold
                //  edges already handled.
                void handleFeatureAngle
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const scalar minCos,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion on warped faces
                void handleWarpedFaces
                (
                    const indirectPrimitivePatch& pp,
                    const scalar faceRatio,
                    const scalar edge0Len,
                    const labelList& cellLevel,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Determine the number of layers per point from the number of
                //  layers per surface.
                void setNumLayers
                (
                    const labelList& patchToNLayers,
                    const labelList& patchIDs,
                    const indirectPrimitivePatch& pp,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Grow no-extrusion layer.
                static void growNoExtrusion
                (
                    const indirectPrimitivePatch& pp,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Calculate pointwise wanted and minimum thickness.
                //  thickness: wanted thickness
                //  minthickness: when to give up and not extrude
                //  Gets per patch parameters and determine pp pointwise
                //  parameters.
                void calculateLayerThickness
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& patchIDs,

                    const scalarField& patchExpansionRatio,
                    const bool relativeSizes,
                    const scalarField& patchFinalLayerThickness,
                    const scalarField& patchMinThickness,

                    const labelList& cellLevel,
                    const labelList& patchNLayers,
                    const scalar edge0Len,

                    scalarField& thickness,
                    scalarField& minThickness,
                    scalarField& expansionRatio
                ) const;


            // Extrusion execution

                //- Synchronize displacement among coupled patches.
                void syncPatchDisplacement
                (
                    const motionSmoother& meshMover,
                    const scalarField& minThickness,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Get nearest point on surface to snap to
                void getPatchDisplacement
                (
                    const motionSmoother& meshMover,
                    const scalarField& thickness,
                    const scalarField& minThickness,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Truncates displacement
                // - for all patchFaces in the faceset displacement gets set
                //   to zero
                // - all displacement < minThickness gets set to zero
                label truncateDisplacement
                (
                    const motionSmoother& meshMover,
                    const scalarField& minThickness,
                    const faceSet& illegalPatchFaces,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Setup layer information (at points and faces) to
                //  modify mesh topology in
                //  regions where layer mesh terminates. Guarantees an
                //  optional slow decreasing of the number of layers.
                //  Returns the number of layers per face and per point
                //  to go into the actual layer addition engine.
                void setupLayerInfoTruncation
                (
                    const motionSmoother& meshMover,
                    const labelList& patchNLayers,
                    const List<extrudeMode>& extrudeStatus,
                    const label nBufferCellsNoExtrude,
                    labelList& nPatchPointLayers,
                    labelList& nPatchFaceLayers
                ) const;

                //- Does any of the cells use a face from faces?
                static bool cellsUseFace
                (
                    const polyMesh& mesh,
                    const labelList& cellLabels,
                    const labelHashSet& faces
                );

                //- Checks the newly added cells and locally unmarks points
                //  so they will not get extruded next time round. Returns
                //  global number of unmarked points (0 if all was fine)
                static label checkAndUnmark
                (
                    const addPatchCellLayer& addLayer,
                    const dictionary& motionDict,
                    const indirectPrimitivePatch& pp,
                    const fvMesh&,

                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Count global number of extruded faces
                static label countExtrusion
                (
                    const indirectPrimitivePatch& pp,
                    const List<extrudeMode>& extrudeStatus
                );

                //- Collect layer faces and layer cells into bools
                //  for ease of handling
                static void getLayerCellsFaces
                (
                    const polyMesh&,
                    const addPatchCellLayer&,
                    boolList&,
                    boolList&
                );

            // Mesh shrinking (to create space for layers)

                //- Average field (over all subset of mesh points) by
                //  summing contribution from edges. Global parallel since only
                //  does master edges for coupled edges.
                template<class Type>
                static void averageNeighbours
                (
                    const polyMesh& mesh,
                    const PackedBoolList& isMasterEdge,
                    const labelList& meshEdges,
                    const labelList& meshPoints,
                    const edgeList& edges,
                    const scalarField& invSumWeight,
                    const Field<Type>& data,
                    Field<Type>& average
                );

                //- Calculate inverse sum of edge weights (currently always 1.0)
                void sumWeights
                (
                    const PackedBoolList& isMasterEdge,
                    const labelList& meshEdges,
                    const labelList& meshPoints,
                    const edgeList& edges,
                    scalarField& invSumWeight
                ) const;

                //- Smooth scalar field on patch
                void smoothField
                (
                    const motionSmoother& meshMover,
                    const PackedBoolList& isMasterEdge,
                    const labelList& meshEdges,
                    const scalarField& fieldMin,
                    const label nSmoothDisp,
                    scalarField& field
                ) const;

                //- Smooth normals on patch.
                void smoothPatchNormals
                (
                    const motionSmoother& meshMover,
                    const PackedBoolList& isMasterEdge,
                    const labelList& meshEdges,
                    const label nSmoothDisp,
                    pointField& normals
                ) const;

                //- Smooth normals in interior.
                void smoothNormals
                (
                    const label nSmoothDisp,
                    const PackedBoolList& isMasterEdge,
                    const labelList& fixedPoints,
                    pointVectorField& normals
                ) const;

                bool isMaxEdge
                (
                    const List<pointData>&,
                    const label edgeI,
                    const scalar minCos
                ) const;

                //- Stop layer growth where mesh wraps around edge with a
                //  large feature angle
                void handleFeatureAngleLayerTerminations
                (
                    const indirectPrimitivePatch& pp,
                    const scalar minCos,
                    List<extrudeMode>& extrudeStatus,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    label& nPointCounter
                ) const;

                //- Find isolated islands (points, edges and faces and
                // layer terminations)
                // in the layer mesh and stop any layer growth at these points.
                void findIsolatedRegions
                (
                    const indirectPrimitivePatch& pp,
                    const PackedBoolList& isMasterEdge,
                    const labelList& meshEdges,
                    const scalar minCosLayerTermination,
                    scalarField& field,
                    List<extrudeMode>& extrudeStatus,
                    pointField& patchDisp,
                    labelList& patchNLayers
                ) const;

                // Calculate medial axis fields
                void medialAxisSmoothingInfo
                (
                    const motionSmoother& meshMover,
                    const label nSmoothNormals,
                    const label nSmoothSurfaceNormals,
                    const scalar minMedianAxisAngleCos,

                    pointVectorField& dispVec,
                    pointScalarField& medialRatio,
                    pointScalarField& medialDist
                ) const;

                //- Main routine to shrink mesh
                void shrinkMeshMedialDistance
                (
                    motionSmoother& meshMover,
                    const dictionary& meshQualityDict,
                    const label nSmoothThickness,
                    const scalar maxThicknessToMedialRatio,
                    const label nAllowableErrors,
                    const label nSnap,
                    const scalar minCosLayerTermination,

                    const scalarField& layerThickness,
                    const scalarField& minThickness,

                    const pointVectorField& dispVec,
                    const pointScalarField& medialRatio,
                    const pointScalarField& medialDist,

                    List<extrudeMode>& extrudeStatus,
                    pointField& patchDisp,
                    labelList& patchNLayers
                ) const;



        //- Disallow default bitwise copy construct
        autoLayerDriver(const autoLayerDriver&);

        //- Disallow default bitwise assignment
        void operator=(const autoLayerDriver&);


public:

    //- Runtime type information
    ClassName("autoLayerDriver");

    // Constructors

        //- Construct from components
        autoLayerDriver(meshRefinement& meshRefiner);


    // Member Functions

            //- Merge patch faces on same cell.
            void mergePatchFacesUndo
            (
                const layerParameters& layerParams,
                const dictionary& motionDict
            );

            //- Add cell layers
            void addLayers
            (
                const layerParameters& layerParams,
                const dictionary& motionDict,
                const labelList& patchIDs,
                const label nAllowableErrors,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor
            );

            //- Add layers according to the dictionary settings
            void doLayers
            (
                const dictionary& shrinkDict,
                const dictionary& motionDict,
                const layerParameters& layerParams,
                const bool preBalance,              // balance before adding?
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "autoLayerDriverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
