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

#include "primitiveMesh.H"
#include "pyramidPointFaceRef.H"
#include "ListOps.H"
#include "mathematicalConstants.H"
#include "SortableList.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::primitiveMesh::closedThreshold_  = 1.0e-6;
Foam::scalar Foam::primitiveMesh::aspectThreshold_  = 1000;
Foam::scalar Foam::primitiveMesh::nonOrthThreshold_ = 70;    // deg
Foam::scalar Foam::primitiveMesh::skewThreshold_    = 4;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::primitiveMesh::checkClosedBoundary(const bool report) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkClosedBoundary("
            << "const bool) const: "
            << "checking whether the boundary is closed" << endl;
    }

    // Loop through all boundary faces and sum up the face area vectors.
    // For a closed boundary, this should be zero in all vector components

    vector sumClosed(vector::zero);
    scalar sumMagClosedBoundary = 0;

    const vectorField& areas = faceAreas();

    for (label faceI = nInternalFaces(); faceI < areas.size(); faceI++)
    {
        sumClosed += areas[faceI];
        sumMagClosedBoundary += mag(areas[faceI]);
    }

    reduce(sumClosed, sumOp<vector>());
    reduce(sumMagClosedBoundary, sumOp<scalar>());

    vector openness = sumClosed/(sumMagClosedBoundary + VSMALL);

    if (cmptMax(cmptMag(openness)) > closedThreshold_)
    {
        if (debug || report)
        {
            Info<< " ***Boundary openness " << openness
                << " possible hole in boundary description."
                << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Boundary openness " << openness << " OK."
                << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkClosedCells
(
    const bool report,
    labelHashSet* setPtr,
    labelHashSet* aspectSetPtr,
    const Vector<label>& meshD
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkClosedCells("
            << "const bool, labelHashSet*, labelHashSet*"
            << ", const Vector<label>&) const: "
            << "checking whether cells are closed" << endl;
    }

    // Check that all cells labels are valid
    const cellList& c = cells();

    label nErrorClosed = 0;

    forAll (c, cI)
    {
        const cell& curCell = c[cI];

        if (min(curCell) < 0 || max(curCell) > nFaces())
        {
            if (setPtr)
            {
                setPtr->insert(cI);
            }

            nErrorClosed++;
        }
    }

    if (nErrorClosed > 0)
    {
        if (debug || report)
        {
            Info<< " ***Cells with invalid face labels found, number of cells "
                << nErrorClosed << endl;
        }

        return true;
    }

    // Loop through cell faces and sum up the face area vectors for each cell.
    // This should be zero in all vector components

    vectorField sumClosed(nCells(), vector::zero);
    vectorField sumMagClosed(nCells(), vector::zero);

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();
    const vectorField& areas = faceAreas();

    forAll (own, faceI)
    {
        // Add to owner
        sumClosed[own[faceI]] += areas[faceI];
        sumMagClosed[own[faceI]] += cmptMag(areas[faceI]);
    }

    forAll (nei, faceI)
    {
        // Subtract from neighbour
        sumClosed[nei[faceI]] -= areas[faceI];
        sumMagClosed[nei[faceI]] += cmptMag(areas[faceI]);
    }

    label nOpen = 0;
    scalar maxOpennessCell = 0;

    label nAspect = 0;
    scalar maxAspectRatio = 0;

    const scalarField& vols = cellVolumes();

    label nDims = 0;
    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (meshD[dir] == 1)
        {
            nDims++;
        }
    }


    // Check the sums
    forAll(sumClosed, cellI)
    {
        scalar maxOpenness = 0;

        for(direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            maxOpenness = max
            (
                maxOpenness,
                mag(sumClosed[cellI][cmpt])
               /(sumMagClosed[cellI][cmpt] + VSMALL)
            );
        }

        maxOpennessCell = max(maxOpennessCell, maxOpenness);

        if (maxOpenness > closedThreshold_)
        {
            if (setPtr)
            {
                setPtr->insert(cellI);
            }

            nOpen++;
        }

        // Calculate the aspect ration as the maximum of Cartesian component
        // aspect ratio to the total area hydraulic area aspect ratio
        scalar minCmpt = VGREAT;
        scalar maxCmpt = -VGREAT;
        for (direction dir = 0; dir < vector::nComponents; dir++)
        {
            if (meshD[dir] == 1)
            {
                minCmpt = min(minCmpt, sumMagClosed[cellI][dir]);
                maxCmpt = max(maxCmpt, sumMagClosed[cellI][dir]);
            }
        }

        scalar aspectRatio = maxCmpt/(minCmpt + VSMALL);
        if (nDims == 3)
        {
            aspectRatio = max
            (
                aspectRatio,
                1.0/6.0*cmptSum(sumMagClosed[cellI])/pow(vols[cellI], 2.0/3.0)
            );
        }

        maxAspectRatio = max(maxAspectRatio, aspectRatio);

        if (aspectRatio > aspectThreshold_)
        {
            if (aspectSetPtr)
            {
                aspectSetPtr->insert(cellI);
            }

            nAspect++;
        }
    }

    reduce(nOpen, sumOp<label>());
    reduce(maxOpennessCell, maxOp<scalar>());

    reduce(nAspect, sumOp<label>());
    reduce(maxAspectRatio, maxOp<scalar>());


    if (nOpen > 0)
    {
        if (debug || report)
        {
            Info<< " ***Open cells found, max cell openness: "
                << maxOpennessCell << ", number of open cells " << nOpen
                << endl;
        }

        return true;
    }

    if (nAspect > 0)
    {
        if (debug || report)
        {
            Info<< " ***High aspect ratio cells found, Max aspect ratio: "
                << maxAspectRatio
                << ", number of cells " << nAspect
                << endl;
        }

        return true;
    }

    if (debug || report)
    {
        Info<< "    Max cell openness = " << maxOpennessCell << " OK." << nl
            << "    Max aspect ratio = " << maxAspectRatio << " OK."
            << endl;
    }

    return false;
}


bool Foam::primitiveMesh::checkFaceAreas
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceAreas("
            << "const bool, labelHashSet*) const: "
            << "checking face area magnitudes" << endl;
    }

    const scalarField magFaceAreas = mag(faceAreas());

    scalar minArea = GREAT;
    scalar maxArea = -GREAT;

    forAll (magFaceAreas, faceI)
    {
        if (magFaceAreas[faceI] < VSMALL)
        {
            if (setPtr)
            {
                setPtr->insert(faceI);
            }
        }

        minArea = min(minArea, magFaceAreas[faceI]);
        maxArea = max(maxArea, magFaceAreas[faceI]);
    }

    reduce(minArea, minOp<scalar>());
    reduce(maxArea, maxOp<scalar>());

    if (minArea < VSMALL)
    {
        if (debug || report)
        {
            Info<< " ***Zero or negative face area detected.  "
                "Minimum area: " << minArea << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Minumum face area = " << minArea
                << ". Maximum face area = " << maxArea
                << ".  Face area magnitudes OK." << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkCellVolumes
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkCellVolumes("
            << "const bool, labelHashSet*) const: "
            << "checking cell volumes" << endl;
    }

    const scalarField& vols = cellVolumes();

    scalar minVolume = GREAT;
    scalar maxVolume = -GREAT;

    label nNegVolCells = 0;

    forAll (vols, cellI)
    {
        if (vols[cellI] < VSMALL)
        {
            if (setPtr)
            {
                setPtr->insert(cellI);
            }

            nNegVolCells++;
        }

        minVolume = min(minVolume, vols[cellI]);
        maxVolume = max(maxVolume, vols[cellI]);
    }

    reduce(minVolume, minOp<scalar>());
    reduce(maxVolume, maxOp<scalar>());
    reduce(nNegVolCells, sumOp<label>());

    if (minVolume < VSMALL)
    {
        if (debug || report)
        {
            Info<< " ***Zero or negative cell volume detected.  "
                << "Minimum negative volume: " << minVolume
                << ", Number of negative volume cells: " << nNegVolCells
                << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Min volume = " << minVolume
                << ". Max volume = " << maxVolume
                << ".  Total volume = " << gSum(vols)
                << ".  Cell volumes OK." << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkFaceOrthogonality
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceOrthogonality("
            << "const bool, labelHashSet*) const: "
            << "checking mesh non-orthogonality" << endl;
    }

    // for all internal faces check theat the d dot S product is positive
    const vectorField& centres = cellCentres();
    const vectorField& areas = faceAreas();

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    // Severe nonorthogonality threshold
    const scalar severeNonorthogonalityThreshold =
        ::cos(nonOrthThreshold_/180.0*mathematicalConstant::pi);

    scalar minDDotS = GREAT;

    scalar sumDDotS = 0;

    label severeNonOrth = 0;

    label errorNonOrth = 0;

    forAll (nei, faceI)
    {
        vector d = centres[nei[faceI]] - centres[own[faceI]];
        const vector& s = areas[faceI];

        scalar dDotS = (d & s)/(mag(d)*mag(s) + VSMALL);

        if (dDotS < severeNonorthogonalityThreshold)
        {
            if (dDotS > SMALL)
            {
                if (setPtr)
                {
                    setPtr->insert(faceI);
                }

                severeNonOrth++;
            }
            else
            {
                if (setPtr)
                {
                    setPtr->insert(faceI);
                }

                errorNonOrth++;
            }
        }

        if (dDotS < minDDotS)
        {
            minDDotS = dDotS;
        }

        sumDDotS += dDotS;
    }

    reduce(minDDotS, minOp<scalar>());
    reduce(sumDDotS, sumOp<scalar>());
    reduce(severeNonOrth, sumOp<label>());
    reduce(errorNonOrth, sumOp<label>());

    if (debug || report)
    {
        label neiSize = nei.size();
        reduce(neiSize, sumOp<label>());

        if (neiSize > 0)
        {
            if (debug || report)
            {
                Info<< "    Mesh non-orthogonality Max: "
                    << ::acos(minDDotS)/mathematicalConstant::pi*180.0
                    << " average: " <<
                    ::acos(sumDDotS/neiSize)/mathematicalConstant::pi*180.0
                    << endl;
            }
        }

        if (severeNonOrth > 0)
        {
            Info<< "   *Number of severely non-orthogonal faces: "
                << severeNonOrth << "." << endl;
        }
    }

    if (errorNonOrth > 0)
    {
        if (debug || report)
        {
            Info<< " ***Number of non-orthogonality errors: "
                << errorNonOrth << "." << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Non-orthogonality check OK." << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkFacePyramids
(
    const bool report,
    const scalar minPyrVol,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFacePyramids("
            << "const bool, const scalar, labelHashSet*) const: "
            << "checking face orientation" << endl;
    }

    // check whether face area vector points to the cell with higher label
    const vectorField& ctrs = cellCentres();

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    const faceList& f = faces();

    const pointField& p = points();

    label nErrorPyrs = 0;

    forAll (f, faceI)
    {
        // Create the owner pyramid - it will have negative volume
        scalar pyrVol = pyramidPointFaceRef(f[faceI], ctrs[own[faceI]]).mag(p);

        if (pyrVol > -minPyrVol)
        {
            if (setPtr)
            {
                setPtr->insert(faceI);
            }

            nErrorPyrs++;
        }

        if (isInternalFace(faceI))
        {
            // Create the neighbour pyramid - it will have positive volume
            scalar pyrVol =
                pyramidPointFaceRef(f[faceI], ctrs[nei[faceI]]).mag(p);

            if (pyrVol < minPyrVol)
            {
                if (setPtr)
                {
                    setPtr->insert(faceI);
                }

                nErrorPyrs++;
            }
        }
    }

    reduce(nErrorPyrs, sumOp<label>());

    if (nErrorPyrs > 0)
    {
        if (debug || report)
        {
            Info<< " ***Error in face pyramids: "
                << nErrorPyrs << " faces are incorrectly oriented."
                << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Face pyramids OK." << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkFaceSkewness
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceSkewnesss("
            << "const bool, labelHashSet*) const: "
            << "checking face skewness" << endl;
    }

    // Warn if the skew correction vector is more than skewWarning times
    // larger than the face area vector

    const pointField& p = points();
    const faceList& fcs = faces();

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();
    const vectorField& cellCtrs = cellCentres();
    const vectorField& faceCtrs = faceCentres();
    const vectorField& fAreas = faceAreas();

    scalar maxSkew = 0;
    label nWarnSkew = 0;

    forAll(nei, faceI)
    {
        vector Cpf = faceCtrs[faceI] - cellCtrs[own[faceI]];
        vector d = cellCtrs[nei[faceI]] - cellCtrs[own[faceI]];

        // Skewness vector
        vector sv =
            Cpf - ((fAreas[faceI] & Cpf)/((fAreas[faceI] & d) + SMALL))*d;
        vector svHat = sv/(mag(sv) + VSMALL);

        // Normalisation distance calculated as the approximate distance
        // from the face centre to the edge of the face in the direction of
        // the skewness
        scalar fd = 0.2*mag(d) + VSMALL;
        const face& f = fcs[faceI];
        forAll(f, pi)
        {
            fd = max(fd, mag(svHat & (p[f[pi]] - faceCtrs[faceI])));
        }

        // Normalised skewness
        scalar skewness = mag(sv)/fd;

        // Check if the skewness vector is greater than the PN vector.
        // This does not cause trouble but is a good indication of a poor mesh.
        if (skewness > skewThreshold_)
        {
            if (setPtr)
            {
                setPtr->insert(faceI);
            }

            nWarnSkew++;
        }

        if(skewness > maxSkew)
        {
            maxSkew = skewness;
        }
    }


    // Boundary faces: consider them to have only skewness error.
    // (i.e. treat as if mirror cell on other side)

    for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
    {
        vector Cpf = faceCtrs[faceI] - cellCtrs[own[faceI]];

        vector normal = fAreas[faceI];
        normal /= mag(normal) + VSMALL;
        vector d = normal*(normal & Cpf);


        // Skewness vector
        vector sv = Cpf - ((fAreas[faceI]&Cpf)/((fAreas[faceI]&d)+VSMALL))*d;
        vector svHat = sv/(mag(sv) + VSMALL);

        // Normalisation distance calculated as the approximate distance
        // from the face centre to the edge of the face in the direction of
        // the skewness
        scalar fd = 0.4*mag(d) + VSMALL;
        const face& f = fcs[faceI];
        forAll(f, pi)
        {
            fd = max(fd, mag(svHat & (p[f[pi]] - faceCtrs[faceI])));
        }

        // Normalised skewness
        scalar skewness = mag(sv)/fd;

        // Check if the skewness vector is greater than the PN vector.
        // This does not cause trouble but is a good indication of a poor mesh.
        if (skewness > skewThreshold_)
        {
            if (setPtr)
            {
                setPtr->insert(faceI);
            }

            nWarnSkew++;
        }

        if(skewness > maxSkew)
        {
            maxSkew = skewness;
        }
    }


    reduce(maxSkew, maxOp<scalar>());
    reduce(nWarnSkew, sumOp<label>());

    if (nWarnSkew > 0)
    {
        if (debug || report)
        {
            Info<< " ***Max skewness = " << maxSkew
                << ", " << nWarnSkew << " highly skew faces detected"
                   " which may impair the quality of the results"
                << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Max skewness = " << maxSkew << " OK." << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkPoints
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkPoints"
            << "(const bool, labelHashSet*) const: "
            << "checking points" << endl;
    }

    label nFaceErrors = 0;
    label nCellErrors = 0;

    const labelListList& pf = pointFaces();

    forAll (pf, pointI)
    {
        if (pf[pointI].empty())
        {
            if (setPtr)
            {
                setPtr->insert(pointI);
            }

            nFaceErrors++;
        }
    }


    forAll (pf, pointI)
    {
        const labelList& pc = pointCells(pointI);

        if (pc.empty())
        {
            if (setPtr)
            {
                setPtr->insert(pointI);
            }

            nCellErrors++;
        }
    }

    reduce(nFaceErrors, sumOp<label>());
    reduce(nCellErrors, sumOp<label>());

    if (nFaceErrors > 0 || nCellErrors > 0)
    {
        if (debug || report)
        {
            Info<< " ***Unused points found in the mesh, "
                   "number unused by faces: " << nFaceErrors
                << " number unused by cells: " << nCellErrors
                << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Point usage OK." << endl;
        }

        return false;
    }
}


// Check convexity of angles in a face. Allow a slight non-convexity.
// E.g. maxDeg = 10 allows for angles < 190 (or 10 degrees concavity)
// (if truly concave and points not visible from face centre the face-pyramid
//  check in checkMesh will fail)
bool Foam::primitiveMesh::checkFaceAngles
(
    const bool report,
    const scalar maxDeg,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceAngles"
            << "(const bool, const scalar, labelHashSet*) const: "
            << "checking face angles" << endl;
    }

    if (maxDeg < -SMALL || maxDeg > 180+SMALL)
    {
        FatalErrorIn
        (
            "primitiveMesh::checkFaceAngles"
            "(const bool, const scalar, labelHashSet*)"
        )   << "maxDeg should be [0..180] but is now " << maxDeg
            << exit(FatalError);
    }

    const scalar maxSin = Foam::sin(maxDeg/180.0*mathematicalConstant::pi);

    const pointField& p = points();
    const faceList& fcs = faces();
    vectorField faceNormals(faceAreas());
    faceNormals /= mag(faceNormals) + VSMALL;

    scalar maxEdgeSin = 0.0;

    label nConcave = 0;

    label errorFaceI = -1;

    forAll(fcs, faceI)
    {
        const face& f = fcs[faceI];

        // Get edge from f[0] to f[size-1];
        vector ePrev(p[f[0]] - p[f[f.size()-1]]);
        scalar magEPrev = mag(ePrev);
        ePrev /= magEPrev + VSMALL;

        forAll(f, fp0)
        {
            // Get vertex after fp
            label fp1 = f.fcIndex(fp0);

            // Normalized vector between two consecutive points
            vector e10(p[f[fp1]] - p[f[fp0]]);
            scalar magE10 = mag(e10);
            e10 /= magE10 + VSMALL;

            if (magEPrev > SMALL && magE10 > SMALL)
            {
                vector edgeNormal = ePrev ^ e10;
                scalar magEdgeNormal = mag(edgeNormal);

                if (magEdgeNormal < maxSin)
                {
                    // Edges (almost) aligned -> face is ok.
                }
                else
                {
                    // Check normal
                    edgeNormal /= magEdgeNormal;

                    if ((edgeNormal & faceNormals[faceI]) < SMALL)
                    {
                        if (faceI != errorFaceI)
                        {
                            // Count only one error per face.
                            errorFaceI = faceI;
                            nConcave++;
                        }

                        if (setPtr)
                        {
                            setPtr->insert(faceI);
                        }

                        maxEdgeSin = max(maxEdgeSin, magEdgeNormal);
                    }
                }
            }

            ePrev = e10;
            magEPrev = magE10;
        }
    }

    reduce(nConcave, sumOp<label>());
    reduce(maxEdgeSin, maxOp<scalar>());

    if (nConcave > 0)
    {
        scalar maxConcaveDegr =
            Foam::asin(Foam::min(1.0, maxEdgeSin))
           *180.0/mathematicalConstant::pi;

        if (debug || report)
        {
            Info<< "   *There are " << nConcave
                << " faces with concave angles between consecutive"
                << " edges. Max concave angle = " << maxConcaveDegr
                << " degrees." << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    All angles in faces OK." << endl;
        }

        return false;
    }
}


// Check warpage of faces. Is calculated as the difference between areas of
// individual triangles and the overall area of the face (which ifself is
// is the average of the areas of the individual triangles).
bool Foam::primitiveMesh::checkFaceFlatness
(
    const bool report,
    const scalar warnFlatness,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceFlatness"
            << "(const bool, const scalar, labelHashSet*) const: "
            << "checking face flatness" << endl;
    }

    if (warnFlatness < 0 || warnFlatness > 1)
    {
        FatalErrorIn
        (
            "primitiveMesh::checkFaceFlatness"
            "(const bool, const scalar, labelHashSet*)"
        )   << "warnFlatness should be [0..1] but is now " << warnFlatness
            << exit(FatalError);
    }


    const pointField& p = points();
    const faceList& fcs = faces();
    const pointField& fctrs = faceCentres();

    // Areas are calculated as the sum of areas. (see
    // primitiveMeshFaceCentresAndAreas.C)
    scalarField magAreas(mag(faceAreas()));

    label nWarped = 0;

    scalar minFlatness = GREAT;
    scalar sumFlatness = 0;
    label nSummed = 0;

    forAll(fcs, faceI)
    {
        const face& f = fcs[faceI];

        if (f.size() > 3 && magAreas[faceI] > VSMALL)
        {
            const point& fc = fctrs[faceI];

            // Calculate the sum of magnitude of areas and compare to magnitude
            // of sum of areas.

            scalar sumA = 0.0;

            forAll(f, fp)
            {
                const point& thisPoint = p[f[fp]];
                const point& nextPoint = p[f.nextLabel(fp)];

                // Triangle around fc.
                vector n = 0.5*((nextPoint - thisPoint)^(fc - thisPoint));
                sumA += mag(n);
            }

            scalar flatness = magAreas[faceI] / (sumA+VSMALL);

            sumFlatness += flatness;
            nSummed++;

            minFlatness = min(minFlatness, flatness);

            if (flatness < warnFlatness)
            {
                nWarped++;

                if (setPtr)
                {
                    setPtr->insert(faceI);
                }
            }
        }
    }


    reduce(nWarped, sumOp<label>());
    reduce(minFlatness, minOp<scalar>());

    reduce(nSummed, sumOp<label>());
    reduce(sumFlatness, sumOp<scalar>());

    if (debug || report)
    {
        if (nSummed > 0)
        {
            Info<< "    Face flatness (1 = flat, 0 = butterfly) : average = "
                << sumFlatness / nSummed << "  min = " << minFlatness << endl;
        }
    }


    if (nWarped> 0)
    {
        if (debug || report)
        {
            Info<< "   *There are " << nWarped
                << " faces with ratio between projected and actual area < "
                << warnFlatness << endl;

            Info<< "    Minimum ratio (minimum flatness, maximum warpage) = "
                << minFlatness << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    All face flatness OK." << endl;
        }

        return false;
    }
}


// Check 1D/2Dness of edges. Gets passed the non-empty directions and
// checks all edges in the mesh whether they:
// - have no component in a non-empty direction or
// - are only in a singe non-empty direction.
// Empty direction info is passed in as a vector of labels (synchronised)
// which are 1 if the direction is non-empty, 0 if it is.
bool Foam::primitiveMesh::checkEdgeAlignment
(
    const bool report,
    const Vector<label>& directions,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkEdgeAlignment("
            << "const bool, const Vector<label>&, labelHashSet*) const: "
            << "checking edge alignment" << endl;
    }

    label nDirs = 0;
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        if (directions[cmpt] == 1)
        {
            nDirs++;
        }
        else if (directions[cmpt] != 0)
        {
            FatalErrorIn
            (
                "primitiveMesh::checkEdgeAlignment"
                "(const bool, const Vector<label>&, labelHashSet*)"
            )   << "directions should contain 0 or 1 but is now " << directions
                << exit(FatalError);
        }
    }

    if (nDirs == vector::nComponents)
    {
        return false;
    }



    const pointField& p = points();
    const faceList& fcs = faces();

    EdgeMap<label> edgesInError;

    forAll(fcs, faceI)
    {
        const face& f = fcs[faceI];

        forAll(f, fp)
        {
            label p0 = f[fp];
            label p1 = f.nextLabel(fp);
            if (p0 < p1)
            {
                vector d(p[p1]-p[p0]);
                scalar magD = mag(d);

                if (magD > ROOTVSMALL)
                {
                    d /= magD;

                    // Check how many empty directions are used by the edge.
                    label nEmptyDirs = 0;
                    label nNonEmptyDirs = 0;
                    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
                    {
                        if (mag(d[cmpt]) > 1e-6)
                        {
                            if (directions[cmpt] == 0)
                            {
                                nEmptyDirs++;
                            }
                            else
                            {
                                nNonEmptyDirs++;
                            }
                        }
                    }

                    if (nEmptyDirs == 0)
                    {
                        // Purely in ok directions.
                    }
                    else if (nEmptyDirs == 1)
                    {
                        // Ok if purely in empty directions.
                        if (nNonEmptyDirs > 0)
                        {
                            edgesInError.insert(edge(p0, p1), faceI);
                        }
                    }
                    else if (nEmptyDirs > 1)
                    {
                        // Always an error
                        edgesInError.insert(edge(p0, p1), faceI);
                    }
                }
            }
        }
    }

    label nErrorEdges = returnReduce(edgesInError.size(), sumOp<label>());

    if (nErrorEdges > 0)
    {
        if (debug || report)
        {
            Info<< " ***Number of edges not aligned with or perpendicular to "
                << "non-empty directions: " << nErrorEdges << endl;
        }

        if (setPtr)
        {
            setPtr->resize(2*edgesInError.size());
            forAllConstIter(EdgeMap<label>, edgesInError, iter)
            {
                setPtr->insert(iter.key()[0]);
                setPtr->insert(iter.key()[1]);
            }
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    All edges aligned with or perpendicular to "
                << "non-empty directions." << endl;
        }
        return false;
    }
}


bool Foam::primitiveMesh::checkUpperTriangular
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkUpperTriangular("
            << "const bool, labelHashSet*) const: "
            << "checking face ordering" << endl;
    }

    // Check whether internal faces are ordered in the upper triangular order
    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    const cellList& c = cells();

    label internal = nInternalFaces();

    // Has error occurred?
    bool error = false;
    // Have multiple faces been detected?
    label nMultipleCells = false;

    // Loop through faceCells once more and make sure that for internal cell
    // the first label is smaller
    for (label faceI = 0; faceI < internal; faceI++)
    {
        if (own[faceI] >= nei[faceI])
        {
            error  = true;

            if (setPtr)
            {
                setPtr->insert(faceI);
            }
        }
    }

    // Loop through all cells. For each cell, find the face that is internal
    // and add it to the check list (upper triangular order).
    // Once the list is completed, check it against the faceCell list

    forAll (c, cellI)
    {
        const labelList& curFaces = c[cellI];

        // Neighbouring cells
        SortableList<label> nbr(curFaces.size());

        forAll(curFaces, i)
        {
            label faceI = curFaces[i];

            if (faceI >= nInternalFaces())
            {
                // Sort last
                nbr[i] = labelMax;
            }
            else
            {
                label nbrCellI = nei[faceI];

                if (nbrCellI == cellI)
                {
                    nbrCellI = own[faceI];
                }

                if (cellI < nbrCellI)
                {
                    // cellI is master
                    nbr[i] = nbrCellI;
                }
                else
                {
                    // nbrCell is master. Let it handle this face.
                    nbr[i] = labelMax;
                }
            }
        }

        nbr.sort();

        // Now nbr holds the cellCells in incremental order. Check:
        // - neighbouring cells appear only once. Since nbr is sorted this
        //   is simple check on consecutive elements
        // - faces indexed in same order as nbr are incrementing as well.

        label prevCell = nbr[0];
        label prevFace = curFaces[nbr.indices()[0]];

        bool hasMultipleFaces = false;

        for (label i = 1; i < nbr.size(); i++)
        {
            label thisCell = nbr[i];
            label thisFace = curFaces[nbr.indices()[i]];

            if (thisCell == labelMax)
            {
                break;
            }

            if (thisCell == prevCell)
            {
                hasMultipleFaces = true;

                if (setPtr)
                {
                    setPtr->insert(prevFace);
                    setPtr->insert(thisFace);
                }
            }
            else if (thisFace < prevFace)
            {
                error = true;

                if (setPtr)
                {
                    setPtr->insert(thisFace);
                }
            }

            prevCell = thisCell;
            prevFace = thisFace;
        }

        if (hasMultipleFaces)
        {
            nMultipleCells++;
        }
    }

    reduce(error, orOp<bool>());
    reduce(nMultipleCells, sumOp<label>());

    if ((debug || report) && nMultipleCells > 0)
    {
        Info<< "  <<Found " << nMultipleCells
            << " neighbouring cells with multiple inbetween faces." << endl;
    }

    if (error)
    {
        if (debug || report)
        {
            Info<< " ***Faces not in upper triangular order." << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Upper triangular ordering OK." << endl;
        }

        return false;
    }
}


bool Foam::primitiveMesh::checkCellsZipUp
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkCellsZipUp("
            << "const bool, labelHashSet*) const: "
            << "checking topological cell openness" << endl;
    }

    label nOpenCells = 0;

    const faceList& f = faces();
    const cellList& c = cells();

    forAll (c, cellI)
    {
        const labelList& curFaces = c[cellI];

        const edgeList cellEdges = c[cellI].edges(f);

        labelList edgeUsage(cellEdges.size(), 0);

        forAll (curFaces, faceI)
        {
            edgeList curFaceEdges = f[curFaces[faceI]].edges();

            forAll (curFaceEdges, faceEdgeI)
            {
                const edge& curEdge = curFaceEdges[faceEdgeI];

                forAll (cellEdges, cellEdgeI)
                {
                    if (cellEdges[cellEdgeI] == curEdge)
                    {
                        edgeUsage[cellEdgeI]++;
                        break;
                    }
                }
            }
        }

        edgeList singleEdges(cellEdges.size());
        label nSingleEdges = 0;

        forAll (edgeUsage, edgeI)
        {
            if (edgeUsage[edgeI] == 1)
            {
                singleEdges[nSingleEdges] = cellEdges[edgeI];
                nSingleEdges++;
            }
            else if (edgeUsage[edgeI] != 2)
            {
                if (setPtr)
                {
                    setPtr->insert(cellI);
                }
            }
        }

        if (nSingleEdges > 0)
        {
            if (setPtr)
            {
                setPtr->insert(cellI);
            }

            nOpenCells++;
        }
    }

    reduce(nOpenCells, sumOp<label>());

    if (nOpenCells > 0)
    {
        if (debug || report)
        {
            Info<< " ***Open cells found, number of cells: " << nOpenCells
                << ". This problem may be fixable using the zipUpMesh utility."
                << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Topological cell zip-up check OK." << endl;
        }

        return false;
    }
}


// Vertices of face within point range and unique.
bool Foam::primitiveMesh::checkFaceVertices
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceVertices("
            << "const bool, labelHashSet*) const: "
            << "checking face vertices" << endl;
    }

    // Check that all vertex labels are valid
    const faceList& f = faces();

    label nErrorFaces = 0;

    forAll (f, fI)
    {
        const face& curFace = f[fI];

        if (min(curFace) < 0 || max(curFace) > nPoints())
        {
            if (setPtr)
            {
                setPtr->insert(fI);
            }

            nErrorFaces++;
        }

        // Uniqueness of vertices
        labelHashSet facePoints(2*curFace.size());

        forAll(curFace, fp)
        {
            bool inserted = facePoints.insert(curFace[fp]);

            if (!inserted)
            {
                if (setPtr)
                {
                    setPtr->insert(fI);
                }

                nErrorFaces++;
            }
        }
    }

    reduce(nErrorFaces, sumOp<label>());

    if (nErrorFaces > 0)
    {
        if (debug || report)
        {
            Info<< "    Faces with invalid vertex labels found, "
                << " number of faces: " << nErrorFaces << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Face vertices OK." << endl;
        }

        return false;
    }
}


// Check if all points on face are shared between faces.
bool Foam::primitiveMesh::checkDuplicateFaces
(
    const label faceI,
    const Map<label>& nCommonPoints,
    label& nBaffleFaces,
    labelHashSet* setPtr
) const
{
    bool error = false;

    forAllConstIter(Map<label>, nCommonPoints, iter)
    {
        label nbFaceI = iter.key();
        label nCommon = iter();

        const face& curFace = faces()[faceI];
        const face& nbFace = faces()[nbFaceI];

        if (nCommon == nbFace.size() || nCommon == curFace.size())
        {
            if (nbFace.size() != curFace.size())
            {
                error = true;
            }
            else
            {
                nBaffleFaces++;
            }

            if (setPtr)
            {
                setPtr->insert(faceI);
                setPtr->insert(nbFaceI);
            }
        }
    }

    return error;
}


// Check that shared points are in consecutive order.
bool Foam::primitiveMesh::checkCommonOrder
(
    const label faceI,
    const Map<label>& nCommonPoints,
    labelHashSet* setPtr
) const
{
    bool error = false;

    forAllConstIter(Map<label>, nCommonPoints, iter)
    {
        label nbFaceI = iter.key();
        label nCommon = iter();

        const face& curFace = faces()[faceI];
        const face& nbFace = faces()[nbFaceI];

        if
        (
            nCommon >= 2
         && nCommon != nbFace.size()
         && nCommon != curFace.size()
        )
        {
            forAll(curFace, fp)
            {
                // Get the index in the neighbouring face shared with curFace
                label nb = findIndex(nbFace, curFace[fp]);

                if (nb != -1)
                {

                    // Check the whole face from nb onwards for shared vertices
                    // with neighbouring face. Rule is that any shared vertices
                    // should be consecutive on both faces i.e. if they are
                    // vertices fp,fp+1,fp+2 on one face they should be
                    // vertices nb, nb+1, nb+2 (or nb+2, nb+1, nb) on the
                    // other face.


                    // Vertices before and after on curFace
                    label fpPlus1 = curFace.fcIndex(fp);
                    label fpMin1  = curFace.rcIndex(fp);

                    // Vertices before and after on nbFace
                    label nbPlus1 = nbFace.fcIndex(nb);
                    label nbMin1  = nbFace.rcIndex(nb);

                    // Find order of walking by comparing next points on both
                    // faces.
                    label curInc = labelMax;
                    label nbInc = labelMax;

                    if (nbFace[nbPlus1] == curFace[fpPlus1])
                    {
                        curInc = 1;
                        nbInc = 1;
                    }
                    else if (nbFace[nbPlus1] == curFace[fpMin1])
                    {
                        curInc = -1;
                        nbInc = 1;
                    }
                    else if (nbFace[nbMin1] == curFace[fpMin1])
                    {
                        curInc = -1;
                        nbInc = -1;
                    }
                    else
                    {
                        curInc = 1;
                        nbInc = -1;
                    }


                    // Pass1: loop until start of common vertices found.
                    label curNb = nb;
                    label curFp = fp;

                    do
                    {
                        curFp += curInc;

                        if (curFp >= curFace.size())
                        {
                            curFp = 0;
                        }
                        else if (curFp < 0)
                        {
                            curFp = curFace.size()-1;
                        }

                        curNb += nbInc;

                        if (curNb >= nbFace.size())
                        {
                            curNb = 0;
                        }
                        else if (curNb < 0)
                        {
                            curNb = nbFace.size()-1;
                        }
                    } while (curFace[curFp] == nbFace[curNb]);


                    // Pass2: check equality walking from curFp, curNb
                    // in opposite order.

                    curInc = -curInc;
                    nbInc = -nbInc;

                    for (label commonI = 0; commonI < nCommon; commonI++)
                    {
                        curFp += curInc;

                        if (curFp >= curFace.size())
                        {
                            curFp = 0;
                        }
                        else if (curFp < 0)
                        {
                            curFp = curFace.size()-1;
                        }

                        curNb += nbInc;

                        if (curNb >= nbFace.size())
                        {
                            curNb = 0;
                        }
                        else if (curNb < 0)
                        {
                            curNb = nbFace.size()-1;
                        }

                        if (curFace[curFp] != nbFace[curNb])
                        {
                            if (setPtr)
                            {
                                setPtr->insert(faceI);
                                setPtr->insert(nbFaceI);
                            }

                            error = true;

                            break;
                        }
                    }


                    // Done the curFace - nbFace combination.
                    break;
                }
            }
        }
    }

    return error;
}


// Checks common vertices between faces. If more than 2 they should be
// consecutive on both faces.
bool Foam::primitiveMesh::checkFaceFaces
(
    const bool report,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkFaceFaces(const bool, labelHashSet*)"
            << " const: " << "checking face-face connectivity" << endl;
    }

    const labelListList& pf = pointFaces();

    label nBaffleFaces = 0;
    label nErrorDuplicate = 0;
    label nErrorOrder = 0;
    Map<label> nCommonPoints(100);

    for (label faceI = 0; faceI < nFaces(); faceI++)
    {
        const face& curFace = faces()[faceI];

        // Calculate number of common points between current faceI and
        // neighbouring face. Store on map.
        nCommonPoints.clear();

        forAll(curFace, fp)
        {
            label pointI = curFace[fp];

            const labelList& nbs = pf[pointI];

            forAll(nbs, nbI)
            {
                label nbFaceI = nbs[nbI];

                if (faceI < nbFaceI)
                {
                    // Only check once for each combination of two faces.

                    Map<label>::iterator fnd = nCommonPoints.find(nbFaceI);

                    if (fnd == nCommonPoints.end())
                    {
                        // First common vertex found.
                        nCommonPoints.insert(nbFaceI, 1);
                    }
                    else
                    {
                        fnd()++;
                    }
                }
            }
        }

        // Perform various checks on common points

        // Check all vertices shared (duplicate point)
        if (checkDuplicateFaces(faceI, nCommonPoints, nBaffleFaces, setPtr))
        {
            nErrorDuplicate++;
        }

        // Check common vertices are consecutive on both faces
        if (checkCommonOrder(faceI, nCommonPoints, setPtr))
        {
            nErrorOrder++;
        }
    }

    reduce(nBaffleFaces, sumOp<label>());
    reduce(nErrorDuplicate, sumOp<label>());
    reduce(nErrorOrder, sumOp<label>());

    if (nBaffleFaces)
    {
        Info<< "    Number of identical duplicate faces (baffle faces): "
            << nBaffleFaces << endl;
    }

    if (nErrorDuplicate > 0 || nErrorOrder > 0)
    {
        if (nErrorDuplicate > 0)
        {
            Info<< " ***Number of duplicate (not baffle) faces found: "
                << nErrorDuplicate << endl;
        }

        if (nErrorOrder > 0)
        {
            Info<< " ***Number of faces with non-consecutive shared points: "
                << nErrorOrder << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Face-face connectivity OK." << endl;
        }

        return false;
    }
}


// Checks cells with 1 or less internal faces. Give numerical problems.
bool Foam::primitiveMesh::checkCellDeterminant
(
    const bool report,    // report,
    labelHashSet* setPtr, // setPtr
    const Vector<label>& meshD
) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkCellDeterminant(const bool"
            << ", labelHashSet*) const: "
            << "checking for under-determined cells" << endl;
    }

    // Determine number of dimensions and (for 2D) missing dimension
    label nDims = 0;
    label twoD = -1;
    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (meshD[dir] == 1)
        {
            nDims++;
        }
        else
        {
            twoD = dir;
        }
    }


    const cellList& c = cells();

    label nErrorCells = 0;

    scalar minDet = GREAT;
    scalar sumDet = 0;
    label nSummed = 0;

    if (nDims == 1)
    {
        minDet = 1;
        sumDet = c.size()*minDet;
        nSummed = c.size();
    }
    else
    {
        forAll (c, cellI)
        {
            const labelList& curFaces = c[cellI];

            // Calculate local normalization factor
            scalar avgArea = 0;

            label nInternalFaces = 0;

            forAll(curFaces, i)
            {
                if (isInternalFace(curFaces[i]))
                {
                    avgArea += mag(faceAreas()[curFaces[i]]);

                    nInternalFaces++;
                }
            }

            if (nInternalFaces == 0)
            {
                if (setPtr)
                {
                    setPtr->insert(cellI);
                }

                nErrorCells++;
            }
            else
            {
                avgArea /= nInternalFaces;

                symmTensor areaTensor(symmTensor::zero);

                forAll(curFaces, i)
                {
                    if (isInternalFace(curFaces[i]))
                    {
                        areaTensor += sqr(faceAreas()[curFaces[i]]/avgArea);
                    }
                }

                if (nDims == 2)
                {
                    // Add the missing eigenvector (such that it does not
                    // affect the determinant)
                    if (twoD == 0)
                    {
                        areaTensor.xx() = 1;
                    }
                    else if (twoD == 1)
                    {
                        areaTensor.yy() = 1;
                    }
                    else
                    {
                        areaTensor.zz() = 1;
                    }
                }

                scalar determinant = mag(det(areaTensor));

                minDet = min(determinant, minDet);
                sumDet += determinant;
                nSummed++;

                if (determinant < 1e-3)
                {
                    if (setPtr)
                    {
                        setPtr->insert(cellI);
                    }

                    nErrorCells++;
                }
            }
        }
    }

    reduce(nErrorCells, sumOp<label>());
    reduce(minDet, minOp<scalar>());
    reduce(sumDet, sumOp<scalar>());
    reduce(nSummed, sumOp<label>());

    if (debug || report)
    {
        if (nSummed > 0)
        {
            Info<< "    Cell determinant (wellposedness) : minimum: " << minDet
                << " average: " << sumDet/nSummed
                << endl;
        }
    }

    if (nErrorCells > 0)
    {
        if (debug || report)
        {
            Info<< " ***Cells with small determinant found, number of cells: "
                << nErrorCells << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Cell determinant check OK." << endl;
        }

        return false;
    }

    return false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::primitiveMesh::checkTopology(const bool report) const
{
    label noFailedChecks = 0;

    if (checkPoints(report)) noFailedChecks++;
    if (checkUpperTriangular(report)) noFailedChecks++;
    if (checkCellsZipUp(report)) noFailedChecks++;
    if (checkFaceVertices(report)) noFailedChecks++;
    if (checkFaceFaces(report)) noFailedChecks++;
    //if (checkCellDeterminant(report)) noFailedChecks++;

    if (noFailedChecks == 0)
    {
        if (debug || report)
        {
            Info<< "    Mesh topology OK." << endl;
        }

        return false;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Failed " << noFailedChecks
                << " mesh topology checks." << endl;
        }

        return true;
    }
}


bool Foam::primitiveMesh::checkGeometry(const bool report) const
{
    label noFailedChecks = 0;

    if (checkClosedBoundary(report)) noFailedChecks++;
    if (checkClosedCells(report)) noFailedChecks++;
    if (checkFaceAreas(report)) noFailedChecks++;
    if (checkCellVolumes(report)) noFailedChecks++;
    if (checkFaceOrthogonality(report)) noFailedChecks++;
    if (checkFacePyramids(report)) noFailedChecks++;
    if (checkFaceSkewness(report)) noFailedChecks++;

    if (noFailedChecks == 0)
    {
        if (debug || report)
        {
            Info<< "    Mesh geometry OK." << endl;
        }
        return false;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Failed " << noFailedChecks
                << " mesh geometry checks." << endl;
        }

        return true;
    }
}


bool Foam::primitiveMesh::checkMesh(const bool report) const
{
    if (debug)
    {
        Info<< "bool primitiveMesh::checkMesh(const bool report) const: "
            << "checking primitiveMesh" << endl;
    }

    label noFailedChecks = checkTopology(report) + checkGeometry(report);

    if (noFailedChecks == 0)
    {
        if (debug || report)
        {
            Info<< "Mesh OK." << endl;
        }

        return false;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Failed " << noFailedChecks
                << " mesh checks." << endl;
        }

        return true;
    }
}


Foam::scalar Foam::primitiveMesh::setClosedThreshold(const scalar val)
{
    scalar prev = closedThreshold_;
    closedThreshold_ = val;

    return prev;
}


Foam::scalar Foam::primitiveMesh::setAspectThreshold(const scalar val)
{
    scalar prev = aspectThreshold_;
    aspectThreshold_ = val;

    return prev;
}


Foam::scalar Foam::primitiveMesh::setNonOrthThreshold(const scalar val)
{
    scalar prev = nonOrthThreshold_;
    nonOrthThreshold_ = val;

    return prev;
}


Foam::scalar Foam::primitiveMesh::setSkewThreshold(const scalar val)
{
    scalar prev = skewThreshold_;
    skewThreshold_ = val;

    return prev;
}


// ************************************************************************* //
