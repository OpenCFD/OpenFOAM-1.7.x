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

Description
    Lightwave OBJ format.

    Note: Java obj loader does not support '#' on line

\*---------------------------------------------------------------------------*/

#include "triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void triSurface::writeOBJ(const bool writeSorted, Ostream& os) const
{
    // Write header
    os  << "# Wavefront OBJ file" << endl
        << "# Regions:" << endl;

    labelList faceMap;

    surfacePatchList myPatches(calcPatches(faceMap));

    const pointField& ps = points();

    // Print patch names as comment
    forAll(myPatches, patchI)
    {
        os  << "#     " << patchI << "    "
            << myPatches[patchI].name() << endl;
    }
    os  << "#" << endl;

    os  << "# points    : " << ps.size() << endl
        << "# triangles : " << size() << endl
        << "#" << endl;


    // Write vertex coords
    forAll(ps, pointi)
    {
        os  << "v "
            << ps[pointi].x() << ' '
            << ps[pointi].y() << ' '
            << ps[pointi].z() << endl;
    }

    if (writeSorted)
    {
        label faceIndex = 0;

        forAll(myPatches, patchI)
        {
            // Print all faces belonging to this patch

            os << "g " << myPatches[patchI].name() << endl;

            for
            (
                label patchFaceI = 0;
                patchFaceI < myPatches[patchI].size();
                patchFaceI++
            )
            {
                const label faceI = faceMap[faceIndex++];

                os  << "f "
                    << operator[](faceI)[0] + 1 << ' '
                    << operator[](faceI)[1] + 1 << ' '
                    << operator[](faceI)[2] + 1
                    //<< "  # " << operator[](faceI).region()
                    << endl;
            }
        }
    }
    else
    {
        forAll(*this, faceI)
        {
            os  << "f "
                << operator[](faceI)[0] + 1 << ' '
                << operator[](faceI)[1] + 1 << ' '
                << operator[](faceI)[2] + 1
                //<< "  # " << operator[](faceI).region()
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
