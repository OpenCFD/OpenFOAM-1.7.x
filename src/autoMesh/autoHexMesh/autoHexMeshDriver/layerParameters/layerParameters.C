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

\*---------------------------------------------------------------------------*/

#include "layerParameters.H"
#include "polyBoundaryMesh.H"
#include "mathematicalConstants.H"
#include "refinementSurfaces.H"
#include "searchableSurfaces.H"
#include "regExp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::layerParameters::defaultConcaveAngle = 90;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Read the number of layers from dictionary. Per patch 0 or the number
// of layers.
Foam::labelList Foam::layerParameters::readNumLayers
(
    const PtrList<dictionary>& surfaceDicts,
    const refinementSurfaces& refineSurfaces,
    const labelList& globalToPatch,
    const polyBoundaryMesh& boundaryMesh
)
{
    // Per surface the number of layers
    labelList globalSurfLayers(surfaceDicts.size());
    // Per surface, per region the number of layers
    List<Map<label> > regionSurfLayers(surfaceDicts.size());

    const labelList& surfaceIndices = refineSurfaces.surfaces();

    forAll(surfaceDicts, surfI)
    {
        const dictionary& dict = surfaceDicts[surfI];

        globalSurfLayers[surfI] = readLabel(dict.lookup("surfaceLayers"));

        if (dict.found("regions"))
        {
            // Per-region layer information

            PtrList<dictionary> regionDicts(dict.lookup("regions"));

            const wordList& regionNames =
                refineSurfaces.geometry()[surfaceIndices[surfI]].regions();

            forAll(regionDicts, dictI)
            {
                const dictionary& regionDict = regionDicts[dictI];

                const word regionName(regionDict.lookup("name"));

                label regionI = findIndex(regionNames, regionName);

                label nLayers = readLabel(regionDict.lookup("surfaceLayers"));

                Info<< "    region " << regionName << ':'<< nl
                    << "        surface layers:" << nLayers << nl;

                regionSurfLayers[surfI].insert(regionI, nLayers);
            }
        }
    }


    // Transfer per surface/region information into patchwise region info

    labelList nLayers(boundaryMesh.size(), 0);

    forAll(surfaceIndices, surfI)
    {
        const wordList& regionNames =
            refineSurfaces.geometry()[surfaceIndices[surfI]].regions();

        forAll(regionNames, regionI)
        {
            const word& regionName = regionNames[regionI];

            label global = refineSurfaces.globalRegion(surfI, regionI);

            label patchI = globalToPatch[global];

            // Initialise to surface-wise layers
            nLayers[patchI] = globalSurfLayers[surfI];

            // Override with region specific data if available
            Map<label>::const_iterator iter =
                regionSurfLayers[surfI].find(regionI);

            if (iter != regionSurfLayers[surfI].end())
            {
                nLayers[patchI] = iter();
            }

            // Check
            if (nLayers[patchI] < 0)
            {
                FatalErrorIn
                (
                    "layerParameters::readNumLayers(..)"
                )   << "Illegal number of layers " << nLayers[patchI]
                    << " for surface "
                    << refineSurfaces.names()[surfI]
                    << " region " << regionName << endl
                    << exit(FatalError);
            }
        }
    }
    return nLayers;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::layerParameters::layerParameters
(
    const PtrList<dictionary>& surfaceDicts,
    const refinementSurfaces& refineSurfaces,
    const labelList& globalToPatch,
    const dictionary& dict,
    const polyBoundaryMesh& boundaryMesh
)
:
    numLayers_
    (
        readNumLayers
        (
            surfaceDicts,
            refineSurfaces,
            globalToPatch,
            boundaryMesh
        )
    ),
    expansionRatio_
    (
        numLayers_.size(),
        readScalar(dict.lookup("expansionRatio"))
    ),
    relativeSizes_(false),
    finalLayerThickness_
    (
        numLayers_.size(),
        readScalar(dict.lookup("finalLayerRatio"))
    ),
    minThickness_
    (
        numLayers_.size(),
        readScalar(dict.lookup("minThickness"))
    ),
    featureAngle_(readScalar(dict.lookup("featureAngle"))),
    concaveAngle_
    (
        dict.lookupOrDefault("concaveAngle", defaultConcaveAngle)
    ),
    nGrow_(readLabel(dict.lookup("nGrow"))),
    nSmoothSurfaceNormals_
    (
        readLabel(dict.lookup("nSmoothSurfaceNormals"))
    ),
    nSmoothNormals_(readLabel(dict.lookup("nSmoothNormals"))),
    nSmoothThickness_(readLabel(dict.lookup("nSmoothThickness"))),
    maxFaceThicknessRatio_
    (
        readScalar(dict.lookup("maxFaceThicknessRatio"))
    ),
    layerTerminationCos_
    (
        Foam::cos
        (
            0.5
          * featureAngle_
          * mathematicalConstant::pi/180.
        )
    ),
    maxThicknessToMedialRatio_
    (
        readScalar(dict.lookup("maxThicknessToMedialRatio"))
    ),
    minMedianAxisAngleCos_
    (
        Foam::cos(readScalar(dict.lookup("minMedianAxisAngle")))
      * mathematicalConstant::pi/180.
    ),
    nBufferCellsNoExtrude_
    (
        readLabel(dict.lookup("nBufferCellsNoExtrude"))
    ),
    nSnap_(readLabel(dict.lookup("nSnap"))),
    nLayerIter_(readLabel(dict.lookup("nLayerIter"))),
    nRelaxedIter_(labelMax)
{
    if (dict.found("nRelaxedIter"))
    {
        dict.lookup("nRelaxedIter") >> nRelaxedIter_;
    }

    if (nLayerIter_ < 0 || nRelaxedIter_ < 0)
    {
        FatalErrorIn("layerParameters::layerParameters(..)")
            << "Layer iterations should be >= 0." << endl
            << "nLayerIter:" << nLayerIter_
            << " nRelaxedIter:" << nRelaxedIter_
            << exit(FatalError);
    }
}


// Construct from dictionary
Foam::layerParameters::layerParameters
(
    const dictionary& dict,
    const polyBoundaryMesh& boundaryMesh
)
:
    numLayers_(boundaryMesh.size(), 0),
    expansionRatio_
    (
        boundaryMesh.size(),
        readScalar(dict.lookup("expansionRatio"))
    ),
    relativeSizes_(dict.lookup("relativeSizes")),
    finalLayerThickness_
    (
        boundaryMesh.size(),
        readScalar(dict.lookup("finalLayerThickness"))
    ),
    minThickness_
    (
        boundaryMesh.size(),
        readScalar(dict.lookup("minThickness"))
    ),
    featureAngle_(readScalar(dict.lookup("featureAngle"))),
    concaveAngle_
    (
        dict.lookupOrDefault("concaveAngle", defaultConcaveAngle)
    ),
    nGrow_(readLabel(dict.lookup("nGrow"))),
    nSmoothSurfaceNormals_
    (
        readLabel(dict.lookup("nSmoothSurfaceNormals"))
    ),
    nSmoothNormals_(readLabel(dict.lookup("nSmoothNormals"))),
    nSmoothThickness_(readLabel(dict.lookup("nSmoothThickness"))),
    maxFaceThicknessRatio_
    (
        readScalar(dict.lookup("maxFaceThicknessRatio"))
    ),
    layerTerminationCos_
    (
        Foam::cos
        (
            0.5
          * featureAngle_
          * mathematicalConstant::pi/180.
        )
    ),
    maxThicknessToMedialRatio_
    (
        readScalar(dict.lookup("maxThicknessToMedialRatio"))
    ),
    minMedianAxisAngleCos_
    (
        Foam::cos(readScalar(dict.lookup("minMedianAxisAngle")))
      * mathematicalConstant::pi/180.
    ),
    nBufferCellsNoExtrude_
    (
        readLabel(dict.lookup("nBufferCellsNoExtrude"))
    ),
    nSnap_(readLabel(dict.lookup("nRelaxIter"))),
    nLayerIter_(readLabel(dict.lookup("nLayerIter"))),
    nRelaxedIter_(labelMax)
{
    if (dict.found("nRelaxedIter"))
    {
        dict.lookup("nRelaxedIter") >> nRelaxedIter_;
    }
    if (nLayerIter_ < 0 || nRelaxedIter_ < 0)
    {
        FatalErrorIn("layerParameters::layerParameters(..)")
            << "Layer iterations should be >= 0." << endl
            << "nLayerIter:" << nLayerIter_
            << " nRelaxedIter:" << nRelaxedIter_
            << exit(FatalError);
    }


    const dictionary& layersDict = dict.subDict("layers");

    forAll(boundaryMesh, patchI)
    {
        const word& patchName = boundaryMesh[patchI].name();

        if (layersDict.found(patchName))
        {
            const dictionary& layerDict = layersDict.subDict(patchName);

            numLayers_[patchI] =
                readLabel(layerDict.lookup("nSurfaceLayers"));

            layerDict.readIfPresent
            (
                "expansionRatio",
                expansionRatio_[patchI]
            );
            layerDict.readIfPresent
            (
                "finalLayerThickness",
                finalLayerThickness_[patchI]
            );
            layerDict.readIfPresent
            (
                "minThickness",
                minThickness_[patchI]
            );
        }
    }


    // Check whether layer specification matches any patches
    const List<keyType> wildCards = layersDict.keys(true);

    forAll(wildCards, i)
    {
        regExp re(wildCards[i]);

        bool hasMatch = false;
        forAll(boundaryMesh, patchI)
        {
            if (re.match(boundaryMesh[patchI].name()))
            {
                hasMatch = true;
                break;
            }
        }
        if (!hasMatch)
        {
            IOWarningIn("layerParameters::layerParameters(..)", layersDict)
                << "Wildcard layer specification for " << wildCards[i]
                << " does not match any patch." << endl;
        }
    }

    const List<keyType> nonWildCards = layersDict.keys(false);

    forAll(nonWildCards, i)
    {
        if (boundaryMesh.findPatchID(nonWildCards[i]) == -1)
        {
            IOWarningIn("layerParameters::layerParameters(..)", layersDict)
                << "Layer specification for " << nonWildCards[i]
                << " does not match any patch." << endl;
        }
    }
}


// ************************************************************************* //
