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

Description
    Given a set of points, find out if the mesh resulting from point motion will
    be valid without actually changing the mesh.

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
#include "pyramidPointFaceRef.H"
#include "cell.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::primitiveMesh::checkMeshMotion
(
    const pointField& newPoints,
    const bool report
) const
{
    if (debug || report)
    {
        Pout<< "bool primitiveMesh::checkMeshMotion("
            << "const pointField& newPoints, const bool report) const: "
            << "checking mesh motion" << endl;
    }

    bool error = false;

    const faceList& f = faces();

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    vectorField fCtrs(nFaces());
    vectorField fAreas(nFaces());

    makeFaceCentresAndAreas(newPoints, fCtrs, fAreas);

    // Check cell volumes and calculate new cell centres
    vectorField cellCtrs(nCells());
    scalarField cellVols(nCells());

    makeCellCentresAndVols(fCtrs, fAreas, cellCtrs, cellVols);

    scalar minVolume = GREAT;
    label nNegVols = 0;

    forAll (cellVols, cellI)
    {
        if (cellVols[cellI] < VSMALL)
        {
            if (debug || report)
            {
                Pout<< "Zero or negative cell volume detected for cell "
                    << cellI << ".  Volume = " << cellVols[cellI] << endl;
            }

            nNegVols++;
        }

        minVolume = min(minVolume, cellVols[cellI]);
    }

    if (nNegVols > 0)
    {
        error = true;

        Pout<< "Zero or negative cell volume in mesh motion in " << nNegVols
            << " cells.  Min volume: " << minVolume << endl;
    }
    else
    {
        if (debug || report)
        {
            Pout<< "Min volume = " << minVolume
                << ".  Total volume = " << sum(cellVols)
                << ".  Cell volumes OK." << endl;
        }
    }

    // Check face areas

    scalar minArea = GREAT;
    label nNegAreas = 0;
    label nPyrErrors = 0;
    label nDotProductErrors = 0;

    forAll (f, faceI)
    {
        const scalar a = Foam::mag(fAreas[faceI]);

        if (a < VSMALL)
        {
            if (debug || report)
            {
                if (isInternalFace(faceI))
                {
                    Pout<< "Zero or negative face area detected for "
                        << "internal face "<< faceI << " between cells "
                        << own[faceI] << " and " << nei[faceI]
                        << ".  Face area magnitude = " << a << endl;
                }
                else
                {
                    Pout<< "Zero or negative face area detected for "
                        << "boundary face " << faceI << " next to cell "
                        << own[faceI] << ".  Face area magnitude = "
                        << a << endl;
                }
            }

            nNegAreas++;
        }

        minArea = min(minArea, a);

        // Create the owner pyramid - it will have negative volume
        scalar pyrVol =
            pyramidPointFaceRef(f[faceI], cellCtrs[own[faceI]]).mag(newPoints);

        if (pyrVol > SMALL)
        {
            if (debug || report)
            {
                Pout<< "Negative pyramid volume: " << -pyrVol
                    << " for face " << faceI << " " << f[faceI]
                    << "  and owner cell: " << own[faceI] << endl
                    << "Owner cell vertex labels: "
                    << cells()[own[faceI]].labels(f)
                    << endl;
            }

            nPyrErrors++;
        }

        if (isInternalFace(faceI))
        {
            // Create the neighbour pyramid - it will have positive volume
            scalar pyrVol =
                pyramidPointFaceRef
                (
                    f[faceI],
                    cellCtrs[nei[faceI]]
                ).mag(newPoints);

            if (pyrVol < -SMALL)
            {
                if (debug || report)
                {
                    Pout<< "Negative pyramid volume: " << pyrVol
                        << " for face " << faceI << " " << f[faceI]
                        << "  and neighbour cell: " << nei[faceI] << nl
                        << "Neighbour cell vertex labels: "
                        << cells()[nei[faceI]].labels(f)
                        << endl;
                }

                nPyrErrors++;
            }

            const vector d = cellCtrs[nei[faceI]] - cellCtrs[own[faceI]];
            const vector& s = fAreas[faceI];
            scalar dDotS = (d & s)/(mag(d)*mag(s) + VSMALL);

            // Only write full message the first time
            if (dDotS < SMALL && nDotProductErrors == 0)
            {
                // Non-orthogonality greater than 90 deg
                WarningIn
                (
                    "primitiveMesh::checkMeshMotion"
                    "(const pointField& newPoints, const bool report) const"
                )   << "Severe non-orthogonality in mesh motion for face "
                    << faceI
                    << " between cells " << own[faceI] << " and " << nei[faceI]
                    << ": Angle = " << ::acos(dDotS)/mathematicalConstant::pi*180.0
                    << " deg." << endl;

                nDotProductErrors++;
            }
        }
    }

    if (nNegAreas > 0)
    {
        error = true;

        WarningIn
        (
            "primitiveMesh::checkMeshMotion"
            "(const pointField& newPoints, const bool report) const"
        )   << "Zero or negative face area in mesh motion in " << nNegAreas
            << " faces.  Min area: " << minArea << endl;
    }
    else
    {
        if (debug || report)
        {
            Pout<< "Min area = " << minArea
                << ".  Face areas OK." << endl;
        }
    }

    if (nPyrErrors > 0)
    {
        Pout<< "Detected " << nPyrErrors
            << " negative pyramid volume in mesh motion" << endl;

        error = true;
    }
    else
    {
        if (debug || report)
        {
            Pout<< "Pyramid volumes OK." << endl;
        }
    }

    if (nDotProductErrors > 0)
    {
        Pout<< "Detected " << nDotProductErrors
            << " in non-orthogonality in mesh motion." << endl;

        error = true;
    }
    else
    {
        if (debug || report)
        {
            Pout<< "Non-orthogonality check OK." << endl;
        }
    }

    if (!error && (debug || report))
    {
        Pout << "Mesh motion check OK." << endl;
    }

    return error;
}


// ************************************************************************* //
