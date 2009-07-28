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

#include "triSurface.H"
#include "STLtriangle.H"
#include "primitivePatch.H"
#include "HashTable.H"
#include "hashSignedLabel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void triSurface::writeSTLASCII(Ostream& os) const
{
    labelList faceMap;

    surfacePatchList myPatches(calcPatches(faceMap));

    label faceIndex = 0;
    forAll(myPatches, patchI)
    {
        // Print all faces belonging to this region
        const surfacePatch& patch = myPatches[patchI];

        os << "solid " << patch.name() << endl;

        for
        (
            label patchFaceI = 0;
            patchFaceI < patch.size();
            patchFaceI++
        )
        {
            const label faceI = faceMap[faceIndex++];

            const vector& n = faceNormals()[faceI];

            os  << "  facet normal "
                << n.x() << ' ' << n.y() << ' ' << n.z() << endl;
            os  << "    outer loop" << endl;

            const labelledTri& f = (*this)[faceI];
            const point& pa = points()[f[0]];
            const point& pb = points()[f[1]];
            const point& pc = points()[f[2]];

            os  << "       vertex "
                << pa.x() << ' ' << pa.y() << ' ' << pa.z() << endl;
            os  << "       vertex "
                << pb.x() << ' ' << pb.y() << ' ' << pb.z() << endl;
            os  << "       vertex "
                << pc.x() << ' ' << pc.y() << ' ' << pc.z() << endl;
            os
                << "    endloop" << endl;
            os
                << "  endfacet" << endl;
        }

        os << "endsolid " << patch.name() << endl;
    }
}


void triSurface::writeSTLBINARY(std::ostream& os) const
{
    // Write the STL header
    string header("Foam binary STL", STLheaderSize);
    os.write(header.c_str(), STLheaderSize);

    label nTris = size();
    os.write(reinterpret_cast<char*>(&nTris), sizeof(unsigned int));

    const vectorField& normals = faceNormals();

    forAll(*this, faceI)
    {
        const labelledTri& f = (*this)[faceI];

        // Convert vector into STL single precision
        STLpoint n(normals[faceI]);
        STLpoint pa(points()[f[0]]);
        STLpoint pb(points()[f[1]]);
        STLpoint pc(points()[f[2]]);

        STLtriangle stlTri(n, pa, pb, pc, f.region());

        stlTri.write(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
