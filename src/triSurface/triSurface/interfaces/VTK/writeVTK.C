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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void triSurface::writeVTK(const bool writeSorted, Ostream& os) const
{
    // Write header
    os  << "# vtk DataFile Version 2.0" << nl
        << "triSurface" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA"
        << nl;

    const pointField& ps = points();

    os  << "POINTS " << ps.size() << " float" << nl;

    // Write vertex coords
    forAll(ps, pointi)
    {
        if (pointi > 0 && (pointi % 10) == 0)
        {
            os  << nl;
        }
        else
        {
            os  << ' ';
        }
        os  << ps[pointi].x() << ' '
            << ps[pointi].y() << ' '
            << ps[pointi].z();
    }
    os  << nl;

    os  << "POLYGONS " << size() << ' ' << 4*size() << nl;

    labelList faceMap;
    surfacePatchList myPatches(calcPatches(faceMap));

    if (writeSorted)
    {
        label faceIndex = 0;

        forAll(myPatches, patchI)
        {
            // Print all faces belonging to this patch

            for
            (
                label patchFaceI = 0;
                patchFaceI < myPatches[patchI].size();
                patchFaceI++
            )
            {
                if (faceIndex > 0 && (faceIndex % 10) == 0)
                {
                    os  << nl;
                }
                else
                {
                    os  << ' ';
                }

                const label faceI = faceMap[faceIndex++];

                os  << "3 "
                    << operator[](faceI)[0] << ' '
                    << operator[](faceI)[1] << ' '
                    << operator[](faceI)[2];
            }
        }
        os  << nl;


        // Print region numbers

        os  << "CELL_DATA " << size() << nl;
        os  << "FIELD attributes 1" << nl;
        os  << "region 1 " << size() << " float" << nl;

        faceIndex = 0;

        forAll(myPatches, patchI)
        {
            for
            (
                label patchFaceI = 0;
                patchFaceI < myPatches[patchI].size();
                patchFaceI++
            )
            {
                if (faceIndex > 0 && (faceIndex % 10) == 0)
                {
                    os  << nl;
                }
                else
                {
                    os  << ' ';
                }

                const label faceI = faceMap[faceIndex++];

                os  << operator[](faceI).region();
            }
        }
        os  << nl;
    }
    else
    {
        forAll(*this, faceI)
        {
            if (faceI > 0 && (faceI % 10) == 0)
            {
                os  << nl;
            }
            else
            {
                os  << ' ';
            }
            os  << "3 "
                << operator[](faceI)[0] << ' '
                << operator[](faceI)[1] << ' '
                << operator[](faceI)[2];
        }
        os  << nl;

        os  << "CELL_DATA " << size() << nl;
        os  << "FIELD attributes 1" << nl;
        os  << "region 1 " << size() << " float" << nl;

        forAll(*this, faceI)
        {
            if (faceI > 0 && (faceI % 10) == 0)
            {
                os  << nl;
            }
            else
            {
                os  << ' ';
            }
            os  << operator[](faceI).region();
        }
        os  << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
