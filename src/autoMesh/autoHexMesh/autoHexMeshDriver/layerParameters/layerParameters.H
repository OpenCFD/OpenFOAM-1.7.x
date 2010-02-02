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
    Foam::layerParameters

Description
    Simple container to keep together layer specific information.

SourceFiles
    layerParameters.C

\*---------------------------------------------------------------------------*/

#ifndef layerParameters_H
#define layerParameters_H

#include "dictionary.H"
#include "scalarField.H"
#include "labelList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyBoundaryMesh;
class refinementSurfaces;

/*---------------------------------------------------------------------------*\
                           Class layerParameters Declaration
\*---------------------------------------------------------------------------*/

class layerParameters
{
    // Static data members

        //- Default angle for faces to be convcave
        static const scalar defaultConcaveAngle;


    // Private data

        // Per patch (not region!) information

            //- How many layers to add.
            labelList numLayers_;

            scalarField expansionRatio_;

            Switch relativeSizes_;

            scalarField finalLayerThickness_;

            scalarField minThickness_;


        scalar featureAngle_;

        scalar concaveAngle_;

        label nGrow_;

        label nSmoothSurfaceNormals_;

        label nSmoothNormals_;

        label nSmoothThickness_;

        scalar maxFaceThicknessRatio_;

        scalar layerTerminationCos_;

        scalar maxThicknessToMedialRatio_;

        scalar minMedianAxisAngleCos_;

        label nBufferCellsNoExtrude_;

        label nSnap_;

        label nLayerIter_;

        label nRelaxedIter_;


    // Private Member Functions

        //- Extract patch-wise number of layers
        static labelList readNumLayers
        (
            const PtrList<dictionary>& surfaceDicts,
            const refinementSurfaces& refineSurfaces,
            const labelList& globalToPatch,
            const polyBoundaryMesh& boundaryMesh
        );

        //- Disallow default bitwise copy construct
        layerParameters(const layerParameters&);

        //- Disallow default bitwise assignment
        void operator=(const layerParameters&);


public:

    // Constructors

        //- Construct from dictionary - old syntax
        layerParameters
        (
            const PtrList<dictionary>& surfaceDicts,
            const refinementSurfaces& refineSurfaces,
            const labelList& globalToPatch,
            const dictionary& dict,
            const polyBoundaryMesh& boundaryMesh
        );

        //- Construct from dictionary - new syntax
        layerParameters(const dictionary& dict, const polyBoundaryMesh&);


    // Member Functions

        // Access

            // Per patch information

                //- How many layers to add.
                const labelList& numLayers() const
                {
                    return numLayers_;
                }

                // Expansion factor for layer mesh
                const scalarField& expansionRatio() const
                {
                    return expansionRatio_;
                }

                //- Are size parameters relative to inner cell size or
                //  absolute distances.
                bool relativeSizes() const
                {
                    return relativeSizes_;
                }

                //- Wanted thickness of final added cell layer. If multiple
                //  layers is the thickness of the layer furthest away
                //  from the wall (i.e. nearest the original mesh)
                //  If relativeSize() this number is relative to undistorted
                //  size of the cell outside layer.
                const scalarField& finalLayerThickness() const
                {
                    return finalLayerThickness_;
                }

                //- Minimum thickness of cell layer. If for any reason layer
                //  cannot be above minThickness do not add layer.
                //  If relativeSize() this number is relative to undistorted
                //  size of the cell outside layer.
                const scalarField& minThickness() const
                {
                    return minThickness_;
                }


            scalar featureAngle() const
            {
                return featureAngle_;
            }

            scalar concaveAngle() const
            {
                return concaveAngle_;
            }

            //- If points get not extruded do nGrow layers of connected faces
            //  that are not grown. Is used to not do layers at all close to
            //  features.
            label nGrow() const
            {
                return nGrow_;
            }

            //- Number of smoothing iterations of surface normals
            label nSmoothSurfaceNormals() const
            {
                return nSmoothSurfaceNormals_;
            }

            //- Number of smoothing iterations of interior mesh movement
            //  direction
            label nSmoothNormals() const
            {
                return nSmoothNormals_;
            }

            //- Stop layer growth on highly warped cells
            scalar maxFaceThicknessRatio() const
            {
                return maxFaceThicknessRatio_;
            }

            scalar layerTerminationCos() const
            {
                return layerTerminationCos_;
            }

            //- Smooth layer thickness over surface patches
            label nSmoothThickness() const
            {
                return nSmoothThickness_;
            }

            //- Reduce layer growth where ratio thickness to medial
            //  distance is large
            scalar maxThicknessToMedialRatio() const
            {
                return maxThicknessToMedialRatio_;
            }

            //- Angle used to pick up medial axis points
            scalar minMedianAxisAngleCos() const
            {
                return minMedianAxisAngleCos_;
            }

            //- Create buffer region for new layer terminations
            label nBufferCellsNoExtrude() const
            {
                return nBufferCellsNoExtrude_;
            }

            label nSnap() const
            {
                return nSnap_;
            }

            // Overall

                //- Number of overall layer addition iterations
                label nLayerIter() const
                {
                    return nLayerIter_;
                }

                //- Number of iterations after which relaxed motion rules
                //  are to be used.
                label nRelaxedIter() const
                {
                    return nRelaxedIter_;
                }



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
