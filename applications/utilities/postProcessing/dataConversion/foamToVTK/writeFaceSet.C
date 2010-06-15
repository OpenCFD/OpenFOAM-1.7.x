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

#include "writeFaceSet.H"
#include "OFstream.H"
#include "writeFuns.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void writeFaceSet
(
    const bool binary,
    const vtkMesh& vMesh,
    const faceSet& set,
    const fileName& fileName
)
{
    const faceList& faces = vMesh.mesh().faces();

    std::ofstream pStream(fileName.c_str());

    pStream
        << "# vtk DataFile Version 2.0" << std::endl
        << set.name() << std::endl;
    if (binary)
    {
        pStream << "BINARY" << std::endl;
    }
    else
    {
        pStream << "ASCII" << std::endl;
    }
    pStream << "DATASET POLYDATA" << std::endl;


    //------------------------------------------------------------------
    //
    // Write topology
    // 
    //------------------------------------------------------------------


    // Construct primitivePatch of faces in faceSet.

    faceList setFaces(set.size());
    labelList setFaceLabels(set.size());
    label setFaceI = 0;

    for
    (
        faceSet::const_iterator iter = set.begin();
        iter != set.end();
        ++iter
    )
    {
        setFaceLabels[setFaceI] = iter.key();
        setFaces[setFaceI] = faces[iter.key()];
        setFaceI++;
    }
    primitiveFacePatch fp(setFaces, vMesh.mesh().points());


    // Write points and faces as polygons

    pStream << "POINTS " << fp.nPoints() << " float" << std::endl;

    DynamicList<floatScalar> ptField(3*fp.nPoints());

    writeFuns::insert(fp.localPoints(), ptField);

    writeFuns::write(pStream, binary, ptField);


    label nFaceVerts = 0;

    forAll(fp.localFaces(), faceI)
    {
        nFaceVerts += fp.localFaces()[faceI].size() + 1;
    }
    pStream << "POLYGONS " << fp.size() << ' ' << nFaceVerts
        << std::endl;


    DynamicList<label> vertLabels(nFaceVerts);

    forAll(fp.localFaces(), faceI)
    {
        const face& f = fp.localFaces()[faceI];

        vertLabels.append(f.size());

        writeFuns::insert(f, vertLabels);
    }
    writeFuns::write(pStream, binary, vertLabels);


    //-----------------------------------------------------------------
    //
    // Write data
    // 
    //-----------------------------------------------------------------

    // Write faceID

    pStream
        << "CELL_DATA " << fp.size() << std::endl
        << "FIELD attributes 1" << std::endl;

    // Cell ids first
    pStream << "faceID 1 " << fp.size() << " int" << std::endl;

    writeFuns::write(pStream, binary, setFaceLabels);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
