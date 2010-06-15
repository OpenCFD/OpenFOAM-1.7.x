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

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "STLtriangle.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "gzstream.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool triSurface::readSTLBINARY(const fileName& STLfileName)
{
    bool compressed = false;

    autoPtr<istream> STLfilePtr
    (
        new ifstream(STLfileName.c_str(), std::ios::binary)
    );

    // If the file is compressed, decompress it before reading.
    if (!STLfilePtr->good() && isFile(STLfileName + ".gz", false))
    {
        compressed = true;
        STLfilePtr.reset(new igzstream((STLfileName + ".gz").c_str()));
    }
    istream& STLfile = STLfilePtr();

    if (!STLfile.good())
    {
        FatalErrorIn("triSurface::readSTLBINARY(const fileName&)")
            << "Cannot read file " << STLfileName
            << " or file " << STLfileName + ".gz"
            << exit(FatalError);
    }

    // Read the STL header
    char header[STLheaderSize];
    STLfile.read(header, STLheaderSize);

    // Check that stream is OK, if not this maybe an ASCII file
    if (!STLfile)
    {
        return false;
    }

    // Read the number of triangles in the STl file
    // (note: read as int so we can check whether >2^31)
    int nTris;
    STLfile.read(reinterpret_cast<char*>(&nTris), sizeof(unsigned int));

    // Check that stream is OK and number of triangles is positive,
    // if not this maybe an ASCII file
    if (!STLfile || nTris < 0)
    {
        return false;
    }

    // Compare the size of the file with that expected from the number of tris
    // If the comparison is not sensible then it maybe an ASCII file
    if (!compressed)
    {
        label dataFileSize = Foam::fileSize(STLfileName) - 80;

        if (nTris < dataFileSize/50 || nTris > dataFileSize/25)
        {
            return false;
        }
    }

    // Everything OK so go ahead and read the triangles.

    // Allocate storage for raw points
    pointField rawPoints(3*nTris);

    // Allocate storage for triangles
    setSize(nTris);

    label rawPointI = 0;

    // Read the triangles
    forAll(*this, i)
    {
        // Read an STL triangle
        STLtriangle stlTri(STLfile);

        // Set the rawPoints to the vertices of the STL triangle
        // and set the point labels of the labelledTri
        rawPoints[rawPointI] = stlTri.a();
        operator[](i)[0] = rawPointI++;

        rawPoints[rawPointI] = stlTri.b();
        operator[](i)[1] = rawPointI++;

        rawPoints[rawPointI] = stlTri.c();
        operator[](i)[2] = rawPointI++;

        operator[](i).region() = stlTri.region();
    }

    //STLfile.close();

    stitchTriangles(rawPoints);

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
