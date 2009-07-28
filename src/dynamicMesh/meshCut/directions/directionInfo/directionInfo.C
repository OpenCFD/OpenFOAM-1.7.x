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


\*---------------------------------------------------------------------------*/

#include "directionInfo.H"
#include "hexMatcher.H"
#include "meshTools.H"
#include "polyMesh.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::directionInfo::equal(const edge& e, const label v0, const label v1)
{
    return
        (e.start() == v0 && e.end() == v1)
     || (e.start() == v1 && e.end() == v0);
}


Foam::point Foam::directionInfo::eMid
(
    const primitiveMesh& mesh,
    const label edgeI
)
{
    const edge& e = mesh.edges()[edgeI];

    return
        0.5
      * (mesh.points()[e.start()] + mesh.points()[e.end()]);
}


// Find edge among edgeLabels that uses v0 and v1
Foam::label Foam::directionInfo::findEdge
(
    const primitiveMesh& mesh,
    const labelList& edgeLabels,
    const label v1,
    const label v0
)
{
    forAll(edgeLabels, edgeLabelI)
    {
        label edgeI = edgeLabels[edgeLabelI];

        const edge& e = mesh.edges()[edgeI];

        if (equal(e, v0, v1))
        {
            return edgeI;
        }
    }

    FatalErrorIn("directionInfo::findEdge")
        << "Cannot find an edge among " << edgeLabels << endl
        << "that uses vertices " << v0
        << " and " << v1
        << abort(FatalError);

    return -1;
}


Foam::label Foam::directionInfo::lowest
(
    const label size,
    const label a,
    const label b
)
{
    // Get next point
    label a1 = (a + 1) % size;

    if (a1 == b)
    {
        return a;
    }
    else
    {
        label b1 = (b + 1) % size;

        if (b1 != a)
        {
            FatalErrorIn("directionInfo::lowest")
                << "Problem : a:" << a << " b:" << b << " size:" << size
                << abort(FatalError);
        }

        return b;
    }
}


// Have edge on hex cell. Find corresponding edge on face. Return -1 if none
// found.
Foam::label Foam::directionInfo::edgeToFaceIndex
(
    const primitiveMesh& mesh,
    const label cellI,
    const label faceI,
    const label edgeI
)
{
    if ((edgeI < 0) || (edgeI >= mesh.nEdges()))
    {
        FatalErrorIn("directionInfo::edgeToFaceIndex")
            << "Illegal edge label:" << edgeI
            << " when projecting cut edge from cell " << cellI
            << " to face " << faceI
            << abort(FatalError);
    }

    const edge& e = mesh.edges()[edgeI];

    const face& f = mesh.faces()[faceI];

    // edgeI is either
    // - in faceI. Convert into index in face.
    // - connected (but not in) to face. Return -1.
    // - in face opposite faceI. Convert into index in face.

    label fpA = findIndex(f, e.start());
    label fpB = findIndex(f, e.end());

    if (fpA != -1)
    {
        if (fpB != -1)
        {
            return lowest(f.size(), fpA, fpB);
        }
        else
        {
            // e.start() in face, e.end() not
            return -1;
        }
    }
    else
    {
        if (fpB != -1)
        {
            // e.end() in face, e.start() not
            return -1;
        }
        else
        {
            // Both not in face.
            // e is on opposite face. Determine corresponding edge on this face:
            // - determine two faces using edge (one is the opposite face, 
            //   one is 'side' face
            // - walk on both these faces to opposite edge
            // - check if this opposite edge is on faceI

            label f0I, f1I;

            meshTools::getEdgeFaces(mesh, cellI, edgeI, f0I, f1I);

            // Walk to opposite edge on face f0
            label edge0I =
                meshTools::walkFace(mesh, f0I, edgeI, e.start(), 2);

            // Check if edge on faceI.

            const edge& e0 = mesh.edges()[edge0I];

            fpA = findIndex(f, e0.start());
            fpB = findIndex(f, e0.end());

            if ((fpA != -1) && (fpB != -1))
            {
                return lowest(f.size(), fpA, fpB);
            }

            // Face0 is doesn't have an edge on faceI (so must be the opposite
            // face) so try face1.

            // Walk to opposite edge on face f1
            label edge1I =
                meshTools::walkFace(mesh, f1I, edgeI, e.start(), 2);

            // Check if edge on faceI.
            const edge& e1 = mesh.edges()[edge1I];

            fpA = findIndex(f, e1.start());
            fpB = findIndex(f, e1.end());

            if ((fpA != -1) && (fpB != -1))
            {
                return lowest(f.size(), fpA, fpB);
            }

            FatalErrorIn("directionInfo::edgeToFaceIndex")
                << "Found connected faces " << mesh.faces()[f0I] << " and "
                << mesh.faces()[f1I] << " sharing edge " << edgeI << endl
                << "But none seems to be connected to face " << faceI
                << " vertices:" << f
                << abort(FatalError);

            return -1;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update this cell with neighbouring face information
bool Foam::directionInfo::updateCell
(
    const polyMesh& mesh,
    const label thisCellI,
    const label neighbourFaceI,
    const directionInfo& neighbourInfo,
    const scalar        // tol
)
{
    if (index_ >= -2)
    {
        // Already determined.
        return false;
    }

    if (hexMatcher().isA(mesh, thisCellI))
    {
        const face& f = mesh.faces()[neighbourFaceI];

        if (neighbourInfo.index() == -2)
        {
            // Geometric information from neighbour
            index_ = -2;
        }
        else if (neighbourInfo.index() == -1)
        {
            // Cut tangential to face. Take any edge connected to face
            // but not used in face.

            // Get first edge on face.
            label edgeI = mesh.faceEdges()[neighbourFaceI][0];

            const edge& e = mesh.edges()[edgeI];

            // Find face connected to face through edgeI and on same cell.
            label faceI = 
                meshTools::otherFace
                (
                    mesh,
                    thisCellI,
                    neighbourFaceI,
                    edgeI
                );

            // Find edge on faceI which is connected to e.start() but not edgeI.
            index_ =
                meshTools::otherEdge
                (
                    mesh,
                    mesh.faceEdges()[faceI],
                    edgeI,
                    e.start()
                );
        }
        else
        {
            // Index is a vertex on the face. Convert to mesh edge.

            // Get mesh edge between f[index_] and f[index_+1]
            label v0 = f[neighbourInfo.index()];
            label v1 = f[(neighbourInfo.index() + 1) % f.size()];

            index_ = findEdge(mesh, mesh.faceEdges()[neighbourFaceI], v0, v1);
        }
    }
    else
    {
        // Not a hex so mark this as geometric.
        index_ = -2;
    }


    n_ = neighbourInfo.n();

    return true;
}    


// Update this face with neighbouring cell information
bool Foam::directionInfo::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const label neighbourCellI,
    const directionInfo& neighbourInfo,
    const scalar    // tol
)
{
    // Handle special cases first

    if (index_ >= -2)
    {
        // Already determined
        return false;
    }

    // Handle normal cases where topological or geometrical info comes from
    // neighbouring cell

    if (neighbourInfo.index() >= 0)
    {
        // Neighbour has topological direction (and hence is hex). Find cut
        // edge on face.
        index_ =
            edgeToFaceIndex
            (
                mesh,
                neighbourCellI,
                thisFaceI,
                neighbourInfo.index()
            );
    }
    else
    {
        // Neighbour has geometric information. Use.
        index_ = -2;
    }


    n_ = neighbourInfo.n();

    return true;
}    


// Merge this with information on same face
bool Foam::directionInfo::updateFace
(
    const polyMesh& mesh,
    const label,    // thisFaceI
    const directionInfo& neighbourInfo,
    const scalar    // tol
)
{
    if (index_ >= -2)
    {
        // Already visited.
        return false;
    }
    else
    {
        index_ = neighbourInfo.index();

        n_ = neighbourInfo.n();

        return true;
    }
}    


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::directionInfo& wDist
)
{
    return os << wDist.index_ << wDist.n_;
}


Foam::Istream& Foam::operator>>(Foam::Istream& is, Foam::directionInfo& wDist)
{
    return is >> wDist.index_ >> wDist.n_;
}

// ************************************************************************* //
