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

#include "cellPointWeight.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::cellPointWeight::debug(debug::debugSwitch("cellPointWeight", 0));
Foam::scalar Foam::cellPointWeight::tol(SMALL);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::cellPointWeight::findTetrahedron
(
    const polyMesh& mesh,
    const vector& position,
    const label cellIndex
)
{
    if (debug)
    {
        Pout<< "\nFoam::cellPointWeight::findTetrahedron" << nl
            << "position = " << position << nl
            << "cellIndex = " << cellIndex << endl;
    }

    // Initialise closest triangle variables
    scalar minUVWClose = VGREAT;
    label pointIClose = 0;
    label faceClose = 0;

    const vector& P0 = mesh.cellCentres()[cellIndex];
    const labelList& cellFaces = mesh.cells()[cellIndex];
    const scalar cellVolume = mesh.cellVolumes()[cellIndex];

    // Find the tet that the point occupies
    forAll(cellFaces, faceI)
    {
        // Decompose each face into triangles, making a tet when
        // augmented by the cell centre
        const labelList& facePoints = mesh.faces()[cellFaces[faceI]];

        label pointI = 1;
        while ((pointI + 1) < facePoints.size())
        {
            // Cartesian co-ordinates of the triangle vertices
            const vector& P1 = mesh.points()[facePoints[0]];
            const vector& P2 = mesh.points()[facePoints[pointI]];
            const vector& P3 = mesh.points()[facePoints[pointI + 1]];

            // Edge vectors
            const vector e1 = P1 - P0;
            const vector e2 = P2 - P0;
            const vector e3 = P3 - P0;

            // Solve for interpolation weighting factors

            // Source term
            const vector rhs = position - P0;

            // Determinant of coefficients matrix
            // Note: if det(A) = 0 the tet is degenerate
            const scalar detA =
                e1.x()*e2.y()*e3.z() + e2.x()*e3.y()*e1.z()
              + e3.x()*e1.y()*e2.z() - e1.x()*e3.y()*e2.z()
              - e2.x()*e1.y()*e3.z() - e3.x()*e2.y()*e1.z();

            if (mag(detA/cellVolume) > tol)
            {
                // Solve using Cramers' rule
                const scalar u =
                (
                   rhs.x()*e2.y()*e3.z() + e2.x()*e3.y()*rhs.z()
                  + e3.x()*rhs.y()*e2.z() - rhs.x()*e3.y()*e2.z()
                  - e2.x()*rhs.y()*e3.z() - e3.x()*e2.y()*rhs.z()
                )/detA;

                const scalar v =
                (
                    e1.x()*rhs.y()*e3.z() + rhs.x()*e3.y()*e1.z()
                  + e3.x()*e1.y()*rhs.z() - e1.x()*e3.y()*rhs.z()
                  - rhs.x()*e1.y()*e3.z() - e3.x()*rhs.y()*e1.z()
                )/detA;

                const scalar w =
                (
                    e1.x()*e2.y()*rhs.z() + e2.x()*rhs.y()*e1.z()
                  + rhs.x()*e1.y()*e2.z() - e1.x()*rhs.y()*e2.z()
                  - e2.x()*e1.y()*rhs.z() - rhs.x()*e2.y()*e1.z()
                )/detA;

                // Check if point is in tet
                // value = 0 indicates position lies on a tet face
                if
                (
                   (u + tol > 0) && (v + tol > 0) && (w + tol > 0)
                && (u + v + w < 1 + tol)
                )
                {
                    faceVertices_[0] = facePoints[0];
                    faceVertices_[1] = facePoints[pointI];
                    faceVertices_[2] = facePoints[pointI + 1];

                    weights_[0] = u;
                    weights_[1] = v;
                    weights_[2] = w;
                    weights_[3] = 1.0 - (u + v + w);

                    return;
                }
                else
                {
                    scalar minU = mag(u);
                    scalar minV = mag(v);
                    scalar minW = mag(w);
                    if (minU > 1.0)
                    {
                        minU -= 1.0;
                    }
                    if (minV > 1.0)
                    {
                        minV -= 1.0;
                    }
                    if (minW > 1.0)
                    {
                        minW -= 1.0;
                    }
                    const scalar minUVW = mag(minU + minV + minW);

                    if (minUVW < minUVWClose)
                    {
                        minUVWClose = minUVW;
                        pointIClose = pointI;
                        faceClose = faceI;
                    }
                }
            }

            pointI++;
        }
    }

    if (debug)
    {
        Pout<< "cellPointWeight::findTetrahedron" << nl
            << "    Tetrahedron search failed; using closest tet values to "
            << "point " << nl << "    cell: " << cellIndex << nl << endl;
    }

    const labelList& facePointsClose = mesh.faces()[cellFaces[faceClose]];
    faceVertices_[0] = facePointsClose[0];
    faceVertices_[1] = facePointsClose[pointIClose];
    faceVertices_[2] = facePointsClose[pointIClose + 1];

    weights_[0] = 0.25;
    weights_[1] = 0.25;
    weights_[2] = 0.25;
    weights_[3] = 0.25;
}


void Foam::cellPointWeight::findTriangle
(
    const polyMesh& mesh,
    const vector& position,
    const label faceIndex
)
{
    if (debug)
    {
        Pout<< "\nbool Foam::cellPointWeight::findTriangle" << nl
            << "position = " << position << nl
            << "faceIndex = " << faceIndex << endl;
    }

    // Initialise closest triangle variables
    scalar minUVClose = VGREAT;
    label pointIClose = 0;

    // Decompose each face into triangles, making a tet when
    // augmented by the cell centre
    const labelList& facePoints = mesh.faces()[faceIndex];

    const scalar faceArea2 = magSqr(mesh.faceAreas()[faceIndex]);

    label pointI = 1;
    while ((pointI + 1) < facePoints.size())
    {
        // Cartesian co-ordinates of the triangle vertices
        const vector& P1 = mesh.points()[facePoints[0]];
        const vector& P2 = mesh.points()[facePoints[pointI]];
        const vector& P3 = mesh.points()[facePoints[pointI + 1]];

        // Direction vectors
        vector v1 = position - P1;
        const vector v2 = P2 - P1;
        const vector v3 = P3 - P1;

        // Plane normal
        vector n = v2 ^ v3;
        n /= mag(n);

        // Remove any offset to plane
        v1 -= (n & v1)*v1;

        // Helper variables
        const scalar d12 = v1 & v2;
        const scalar d13 = v1 & v3;
        const scalar d22 = v2 & v2;
        const scalar d23 = v2 & v3;
        const scalar d33 = v3 & v3;

        // Determinant of coefficients matrix
        // Note: if det(A) = 0 the triangle is degenerate
        const scalar detA = d22*d33 - d23*d23;

        if (0.25*detA/faceArea2 > tol)
        {
            // Solve using Cramers' rule
            const scalar u = (d12*d33 - d23*d13)/detA;
            const scalar v = (d22*d13 - d12*d23)/detA;

            // Check if point is in triangle
            if ((u + tol > 0) && (v + tol > 0) && (u + v < 1 + tol))
            {
                // Indices of the cell vertices making up the triangle
                faceVertices_[0] = facePoints[0];
                faceVertices_[1] = facePoints[pointI];
                faceVertices_[2] = facePoints[pointI + 1];

                weights_[0] = u;
                weights_[1] = v;
                weights_[2] = 1.0 - (u + v);
                weights_[3] = 0.0;

                return;
            }
            else
            {
                scalar minU = mag(u);
                scalar minV = mag(v);
                if (minU > 1.0)
                {
                    minU -= 1.0;
                }
                if (minV > 1.0)
                {
                    minV -= 1.0;
                }
                const scalar minUV = mag(minU + minV);

                if (minUV < minUVClose)
                {
                    minUVClose = minUV;
                    pointIClose = pointI;
                }
            }
        }

        pointI++;
    }

    if (debug)
    {
        Pout<< "Foam::cellPointWeight::findTriangle"
            << "Triangle search failed; using closest triangle to point" << nl
            << "    cell face: " << faceIndex << nl << endl;
    }

    // Indices of the cell vertices making up the triangle
    faceVertices_[0] = facePoints[0];
    faceVertices_[1] = facePoints[pointIClose];
    faceVertices_[2] = facePoints[pointIClose + 1];

    weights_[0] = 1.0/3.0;
    weights_[1] = 1.0/3.0;
    weights_[2] = 1.0/3.0;
    weights_[3] = 0.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellPointWeight::cellPointWeight
(
    const polyMesh& mesh,
    const vector& position,
    const label cellIndex,
    const label faceIndex
)
:
    cellIndex_(cellIndex)
{
    if (faceIndex < 0)
    {
        // Face data not supplied
        findTetrahedron(mesh, position, cellIndex);
    }
    else
    {
        // Face data supplied
        findTriangle(mesh, position, faceIndex);
    }
}


// ************************************************************************* //
