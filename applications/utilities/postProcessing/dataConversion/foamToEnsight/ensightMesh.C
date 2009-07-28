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

#include "argList.H"
#include "Time.H"
#include "ensightMesh.H"
#include "fvMesh.H"
#include "globalMeshData.H"
#include "PstreamCombineReduceOps.H"
#include "processorPolyPatch.H"
#include "cellModeller.H"
#include "IOmanip.H"
#include "itoa.H"
#include "ensightWriteBinary.H"
#include <fstream>

// * * * * * * * * * * * * * Private Functions * * * * * * * * * * * * * * //

namespace Foam
{
    //- Proxy-class to hold the patch processor list combination operator
    class concatPatchProcs
    {

    public:

        void operator()
        (
            List<labelList>& x,
            const List<labelList>& y
        ) const
        {
            forAll(y, i)
            {
                const labelList& yPatches = y[i];

                if (yPatches.size())
                {
                    labelList& xPatches = x[i];

                    label offset = xPatches.size();
                    xPatches.setSize(offset + yPatches.size());

                    forAll(yPatches, i)
                    {
                        xPatches[i + offset] = yPatches[i];
                    }
                }
            }
        }
    };
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightMesh::ensightMesh
(
    const fvMesh& mesh,
    const argList& args,
    const bool binary
)
:
    mesh_(mesh),
    binary_(binary),
    patchPartOffset_(2),
    meshCellSets_(mesh_.nCells()),
    boundaryFaceSets_(mesh_.boundary().size()),
    allPatchNames_(0),
    allPatchProcs_(0),
    patchNames_(0),
    nPatchPrims_(0)
{
    const cellShapeList& cellShapes = mesh.cellShapes();

    const cellModel& tet = *(cellModeller::lookup("tet"));
    const cellModel& pyr = *(cellModeller::lookup("pyr"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& wedge = *(cellModeller::lookup("wedge"));
    const cellModel& hex = *(cellModeller::lookup("hex"));

    if (!args.optionFound("noPatches"))
    {
        allPatchNames_ = wordList::subList
        (
            mesh_.boundaryMesh().names(), mesh_.boundary().size()
          - mesh_.globalData().processorPatches().size()
        );

        allPatchProcs_.setSize(allPatchNames_.size());

        forAll (allPatchProcs_, patchi)
        {
            if (mesh_.boundary()[patchi].size())
            {
                allPatchProcs_[patchi].setSize(1);
                allPatchProcs_[patchi][0] = Pstream::myProcNo();
            }
        }

        combineReduce(allPatchProcs_, concatPatchProcs());

        if (args.optionFound("patches"))
        {
            wordList patchNameList(args.optionLookup("patches")());

            if (patchNameList.empty())
            {
                patchNameList = allPatchNames_;
            }

            forAll (patchNameList, i)
            {
                patchNames_.insert(patchNameList[i]);
            }
        }
    }

    if (patchNames_.size())
    {
        // no internalMesh
        patchPartOffset_ = 1;
    }
    else
    {
        // Count the shapes
        labelList& tets = meshCellSets_.tets;
        labelList& pyrs = meshCellSets_.pyrs;
        labelList& prisms = meshCellSets_.prisms;
        labelList& wedges = meshCellSets_.wedges;
        labelList& hexes = meshCellSets_.hexes;
        labelList& polys = meshCellSets_.polys;

        label nTets = 0;
        label nPyrs = 0;
        label nPrisms = 0;
        label nWedges = 0;
        label nHexes = 0;
        label nPolys = 0;

        forAll(cellShapes, cellI)
        {
            const cellShape& cellShape = cellShapes[cellI];
            const cellModel& cellModel = cellShape.model();

            if (cellModel == tet)
            {
                tets[nTets++] = cellI;
            }
            else if (cellModel == pyr)
            {
                pyrs[nPyrs++] = cellI;
            }
            else if (cellModel == prism)
            {
                prisms[nPrisms++] = cellI;
            }
            else if (cellModel == wedge)
            {
                wedges[nWedges++] = cellI;
            }
            else if (cellModel == hex)
            {
                hexes[nHexes++] = cellI;
            }
            else
            {
                polys[nPolys++] = cellI;
            }
        }

        tets.setSize(nTets);
        pyrs.setSize(nPyrs);
        prisms.setSize(nPrisms);
        wedges.setSize(nWedges);
        hexes.setSize(nHexes);
        polys.setSize(nPolys);

        meshCellSets_.nTets = nTets;
        reduce(meshCellSets_.nTets, sumOp<label>());

        meshCellSets_.nPyrs = nPyrs;
        reduce(meshCellSets_.nPyrs, sumOp<label>());

        meshCellSets_.nPrisms = nPrisms;
        reduce(meshCellSets_.nPrisms, sumOp<label>());

        meshCellSets_.nHexesWedges = nHexes + nWedges;
        reduce(meshCellSets_.nHexesWedges, sumOp<label>());

        meshCellSets_.nPolys = nPolys;
        reduce(meshCellSets_.nPolys, sumOp<label>());
    }

    if (!args.optionFound("noPatches"))
    {
        forAll (mesh.boundary(), patchi)
        {
            if (mesh.boundary()[patchi].size())
            {
                const polyPatch& p = mesh.boundaryMesh()[patchi];

                labelList& tris = boundaryFaceSets_[patchi].tris;
                labelList& quads = boundaryFaceSets_[patchi].quads;
                labelList& polys = boundaryFaceSets_[patchi].polys;

                tris.setSize(p.size());
                quads.setSize(p.size());
                polys.setSize(p.size());

                label nTris = 0;
                label nQuads = 0;
                label nPolys = 0;

                forAll(p, faceI)
                {
                    const face& f = p[faceI];

                    if (f.size() == 3)
                    {
                        tris[nTris++] = faceI;
                    }
                    else if (f.size() == 4)
                    {
                        quads[nQuads++] = faceI;
                    }
                    else
                    {
                        polys[nPolys++] = faceI;
                    }
                }

                tris.setSize(nTris);
                quads.setSize(nQuads);
                polys.setSize(nPolys);
            }
        }
    }


    forAll(allPatchNames_, patchi)
    {
        const word& patchName = allPatchNames_[patchi];
        nFacePrimitives nfp;

        if (patchNames_.empty() || patchNames_.found(patchName))
        {
            if (mesh.boundary()[patchi].size())
            {
                nfp.nPoints = mesh.boundaryMesh()[patchi].localPoints().size();
                nfp.nTris   = boundaryFaceSets_[patchi].tris.size();
                nfp.nQuads  = boundaryFaceSets_[patchi].quads.size();
                nfp.nPolys  = boundaryFaceSets_[patchi].polys.size();
            }
        }

        reduce(nfp.nPoints, sumOp<label>());
        reduce(nfp.nTris, sumOp<label>());
        reduce(nfp.nQuads, sumOp<label>());
        reduce(nfp.nPolys, sumOp<label>());

        nPatchPrims_.insert(patchName, nfp);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightMesh::~ensightMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightMesh::writePoints
(
    const scalarField& pointsComponent,
    OFstream& ensightGeometryFile
) const
{
    forAll(pointsComponent, pointI)
    {
        ensightGeometryFile<< setw(12) << float(pointsComponent[pointI]) << nl;
    }
}


Foam::cellShapeList Foam::ensightMesh::map
(
    const cellShapeList& cellShapes,
    const labelList& prims
) const
{
    cellShapeList mcsl(prims.size());

    forAll(prims, i)
    {
        mcsl[i] = cellShapes[prims[i]];
    }

    return mcsl;
}


Foam::cellShapeList Foam::ensightMesh::map
(
    const cellShapeList& cellShapes,
    const labelList& hexes,
    const labelList& wedges
) const
{
    cellShapeList mcsl(hexes.size() + wedges.size());

    forAll(hexes, i)
    {
        mcsl[i] = cellShapes[hexes[i]];
    }

    label offset = hexes.size();

    const cellModel& hex = *(cellModeller::lookup("hex"));
    labelList hexLabels(8);

    forAll(wedges, i)
    {
        const cellShape& cellPoints = cellShapes[wedges[i]];

        hexLabels[0] = cellPoints[0];
        hexLabels[1] = cellPoints[1];
        hexLabels[2] = cellPoints[0];
        hexLabels[3] = cellPoints[2];
        hexLabels[4] = cellPoints[3];
        hexLabels[5] = cellPoints[4];
        hexLabels[6] = cellPoints[6];
        hexLabels[7] = cellPoints[5];

        mcsl[i + offset] = cellShape(hex, hexLabels);
    }

    return mcsl;
}


void Foam::ensightMesh::writePrims
(
    const cellShapeList& cellShapes,
    const label pointOffset,
    OFstream& ensightGeometryFile
) const
{
    label po = pointOffset + 1;

    forAll(cellShapes, i)
    {
        const cellShape& cellPoints = cellShapes[i];

        forAll(cellPoints, pointI)
        {
            ensightGeometryFile<< setw(10) << cellPoints[pointI] + po;
        }
        ensightGeometryFile << nl;
    }
}


void Foam::ensightMesh::writePrimsBinary
(
    const cellShapeList& cellShapes,
    const label pointOffset,
    std::ofstream& ensightGeometryFile
) const
{
    label po = pointOffset + 1;

    // Create a temp int array
    int numElem;

    numElem = cellShapes.size();

    if (cellShapes.size())
    {
        // All the cellShapes have the same number of elements!
        int numIntElem = cellShapes.size()*cellShapes[0].size();
        List<int> temp(numIntElem);

        int n = 0;

        forAll(cellShapes, i)
        {
            const cellShape& cellPoints = cellShapes[i];

            forAll(cellPoints, pointI)
            {
                temp[n] = cellPoints[pointI] + po;
                n++;
            }
        }

        ensightGeometryFile.write
        (
            reinterpret_cast<char*>(temp.begin()),
            numIntElem*sizeof(int)
        );
    }
}


void Foam::ensightMesh::writePolys
(
    const labelList& polys,
    const cellList& cellFaces,
    const faceList& faces,
    const label pointOffset,
    OFstream& ensightGeometryFile
) const
{
    if (polys.size())
    {
        ensightGeometryFile
            << "nfaced" << nl << setw(10) << polys.size() << nl;

        label po = pointOffset + 1;

        forAll(polys, i)
        {
            ensightGeometryFile
                << setw(10) << cellFaces[polys[i]].size() << nl;
        }

        forAll(polys, i)
        {
            const labelList& cf = cellFaces[polys[i]];

            forAll(cf, faceI)
            {
                ensightGeometryFile
                    << setw(10) << faces[cf[faceI]].size() << nl;
            }
        }

        forAll(polys, i)
        {
            const labelList& cf = cellFaces[polys[i]];

            forAll(cf, faceI)
            {
                const face& f = faces[cf[faceI]];

                forAll(f, pointI)
                {
                    ensightGeometryFile << setw(10) << f[pointI] + po;
                }
                ensightGeometryFile << nl;
            }
        }
    }
}


void Foam::ensightMesh::writePolysBinary
(
    const labelList& polys,
    const cellList& cellFaces,
    const faceList& faces,
    const label pointOffset,
    std::ofstream& ensightGeometryFile
) const
{
    if (polys.size())
    {
        writeEnsDataBinary("nfaced",ensightGeometryFile);
        writeEnsDataBinary(polys.size(),ensightGeometryFile);

        label po = pointOffset + 1;

        //TODO No buffer at the moment. To be done for speed purposes!
        forAll(polys, i)
        {
            writeEnsDataBinary
            (
                cellFaces[polys[i]].size(),
                ensightGeometryFile
            );
        }

        forAll(polys, i)
        {
            const labelList& cf = cellFaces[polys[i]];

            forAll(cf, faceI)
            {
                writeEnsDataBinary
                (
                    faces[cf[faceI]].size(),
                    ensightGeometryFile
                );
            }
        }

        forAll(polys, i)
        {
            const labelList& cf = cellFaces[polys[i]];

            forAll(cf, faceI)
            {
                const face& f = faces[cf[faceI]];

                forAll(f, pointI)
                {
                    writeEnsDataBinary(f[pointI] + po,ensightGeometryFile);
                }
            }
        }
    }
}


void Foam::ensightMesh::writeAllPrims
(
    const char* key,
    const label nPrims,
    const cellShapeList& cellShapes,
    const labelList& pointOffsets,
    OFstream& ensightGeometryFile
) const
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            ensightGeometryFile << key << nl << setw(10) << nPrims << nl;

            writePrims(cellShapes, 0, ensightGeometryFile);

            for (int slave=1; slave<Pstream::nProcs(); slave++)
            {
                IPstream fromSlave(Pstream::scheduled, slave);
                cellShapeList cellShapes(fromSlave);

                writePrims
                (
                    cellShapes,
                    pointOffsets[slave-1],
                    ensightGeometryFile
                );
            }
        }
        else
        {
            OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
            toMaster<< cellShapes;
        }
    }
}


void Foam::ensightMesh::writeAllPrimsBinary
(
    const char* key,
    const label nPrims,
    const cellShapeList& cellShapes,
    const labelList& pointOffsets,
    std::ofstream& ensightGeometryFile
) const
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            writeEnsDataBinary(key,ensightGeometryFile);
            writeEnsDataBinary(nPrims,ensightGeometryFile);

            writePrimsBinary(cellShapes, 0, ensightGeometryFile);

            for (int slave=1; slave<Pstream::nProcs(); slave++)
            {
                IPstream fromSlave(Pstream::scheduled, slave);
                cellShapeList cellShapes(fromSlave);

                writePrimsBinary
                (
                    cellShapes,
                    pointOffsets[slave-1],
                    ensightGeometryFile
                );
            }
        }
        else
        {
            OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
            toMaster<< cellShapes;
        }
    }
}


void Foam::ensightMesh::writeFacePrims
(
    const char* key,
    const faceList& patchFaces,
    const label pointOffset,
    OFstream& ensightGeometryFile
) const
{
    if (patchFaces.size())
    {
        if (word(key) == "nsided")
        {
            ensightGeometryFile
                << key << nl << setw(10) << patchFaces.size() << nl;

            forAll(patchFaces, i)
            {
                ensightGeometryFile
                    << setw(10) << patchFaces[i].size() << nl;
            }
        }

        label po = pointOffset + 1;

        forAll(patchFaces, i)
        {
            const face& patchFace = patchFaces[i];

            forAll(patchFace, pointI)
            {
                ensightGeometryFile << setw(10) << patchFace[pointI] + po;
            }
            ensightGeometryFile << nl;
        }
    }
}


void Foam::ensightMesh::writeFacePrimsBinary
(
    const char* key,
    const faceList& patchFaces,
    const label pointOffset,
    std::ofstream& ensightGeometryFile
) const
{
    if (patchFaces.size())
    {
        //TODO No buffer at the moment. To be done for speed purposes!
        if (word(key) == "nsided")
        {
            writeEnsDataBinary(key,ensightGeometryFile);
            writeEnsDataBinary(patchFaces.size(),ensightGeometryFile);

            forAll(patchFaces, i)
            {
                writeEnsDataBinary
                (
                    patchFaces[i].size(),
                    ensightGeometryFile
                );
            }
        }

        label po = pointOffset + 1;

        forAll(patchFaces, i)
        {
            const face& patchFace = patchFaces[i];

            forAll(patchFace, pointI)
            {
                writeEnsDataBinary
                (
                    patchFace[pointI] + po,
                    ensightGeometryFile
                );
            }
        }
    }
}


Foam::faceList Foam::ensightMesh::map
(
    const faceList& patchFaces,
    const labelList& prims
) const
{
    faceList ppf(prims.size());

    forAll (prims, i)
    {
        ppf[i] = patchFaces[prims[i]];
    }

    return ppf;
}


void Foam::ensightMesh::writeAllFacePrims
(
    const char* key,
    const labelList& prims,
    const label nPrims,
    const faceList& patchFaces,
    const labelList& pointOffsets,
    const labelList& patchProcessors,
    OFstream& ensightGeometryFile
) const
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            if (word(key) != "nsided")
            {
                ensightGeometryFile << key << nl << setw(10) << nPrims << nl;
            }

            if (&prims != NULL)
            {
                writeFacePrims
                (
                    key,
                    map(patchFaces, prims),
                    0,
                    ensightGeometryFile
                );
            }

            forAll (patchProcessors, i)
            {
                if (patchProcessors[i] != 0)
                {
                    label slave = patchProcessors[i];
                    IPstream fromSlave(Pstream::scheduled, slave);
                    faceList patchFaces(fromSlave);

                    writeFacePrims
                    (
                        key,
                        patchFaces,
                        pointOffsets[i],
                        ensightGeometryFile
                    );
                }
            }
        }
        else if (&prims != NULL)
        {
            OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
            toMaster<< map(patchFaces, prims);
        }
    }
}


void Foam::ensightMesh::writeAllFacePrimsBinary
(
    const char* key,
    const labelList& prims,
    const label nPrims,
    const faceList& patchFaces,
    const labelList& pointOffsets,
    const labelList& patchProcessors,
    std::ofstream& ensightGeometryFile
) const
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            if (word(key) != "nsided")
            {
                writeEnsDataBinary(key,ensightGeometryFile);
                writeEnsDataBinary(nPrims,ensightGeometryFile);
            }

            if (&prims != NULL)
            {
                writeFacePrimsBinary
                (
                    key,
                    map(patchFaces, prims),
                    0,
                    ensightGeometryFile
                );
            }

            forAll (patchProcessors, i)
            {
                if (patchProcessors[i] != 0)
                {
                    label slave = patchProcessors[i];
                    IPstream fromSlave(Pstream::scheduled, slave);
                    faceList patchFaces(fromSlave);

                    writeFacePrimsBinary
                    (
                        key,
                        patchFaces,
                        pointOffsets[i],
                        ensightGeometryFile
                    );
                }
            }
        }
        else if (&prims != NULL)
        {
            OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
            toMaster<< map(patchFaces, prims);
        }
    }
}


void Foam::ensightMesh::write
(
    const fileName& postProcPath,
    const word& prepend,
    const label timeIndex,
    Ostream& ensightCaseFile
) const
{
    if (binary_)
    {
        writeBinary(postProcPath, prepend, timeIndex, ensightCaseFile);
    }
    else
    {
        writeAscii(postProcPath, prepend, timeIndex, ensightCaseFile);
    }
}


void Foam::ensightMesh::writeAscii
(
    const fileName& postProcPath,
    const word& prepend,
    const label timeIndex,
    Ostream& ensightCaseFile
) const
{
    const Time& runTime = mesh_.time();
    const pointField& points = mesh_.points();
    const cellList& cellFaces = mesh_.cells();
    const faceList& faces = mesh_.faces();
    const cellShapeList& cellShapes = mesh_.cellShapes();

    word timeFile = prepend;

    if (timeIndex == 0)
    {
        timeFile += "000.";
    }
    else if (mesh_.moving())
    {
        timeFile += itoa(timeIndex) + '.';
    }

    // set the filename of the ensight file
    fileName ensightGeometryFileName = timeFile + "mesh";

    OFstream *ensightGeometryFilePtr = NULL;
    if (Pstream::master())
    {
        ensightGeometryFilePtr = new OFstream
        (
            postProcPath/ensightGeometryFileName,
            runTime.writeFormat(),
            runTime.writeVersion(),
            IOstream::UNCOMPRESSED
        );
    }

    OFstream& ensightGeometryFile = *ensightGeometryFilePtr;

    if (Pstream::master())
    {
        // Set Format
        ensightGeometryFile.setf
        (
            ios_base::scientific,
            ios_base::floatfield
        );
        ensightGeometryFile.precision(5);

        ensightGeometryFile
            << "EnSight Geometry File" << nl
            << "written from OpenFOAM-" << Foam::FOAMversion << nl
            << "node id assign" << nl
            << "element id assign" << nl;
    }

    labelList pointOffsets(Pstream::nProcs(), 0);

    if (patchNames_.empty())
    {
        label nPoints = points.size();
        Pstream::gather(nPoints, sumOp<label>());

        if (Pstream::master())
        {
            ensightGeometryFile
                << "part" << nl
                << setw(10) << 1 << nl
                << "internalMesh" << nl
                << "coordinates" << nl
                << setw(10) << nPoints
                << endl;

            for (direction d=0; d<vector::nComponents; d++)
            {
                writePoints(points.component(d), ensightGeometryFile);
                pointOffsets[0] = points.size();

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    IPstream fromSlave(Pstream::scheduled, slave);
                    scalarField pointsComponent(fromSlave);
                    writePoints(pointsComponent, ensightGeometryFile);
                    pointOffsets[slave] =
                        pointOffsets[slave-1]
                      + pointsComponent.size();
                }
            }
        }
        else
        {
            for (direction d=0; d<vector::nComponents; d++)
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< points.component(d);
            }
        }

        writeAllPrims
        (
            "hexa8",
            meshCellSets_.nHexesWedges,
            map(cellShapes, meshCellSets_.hexes, meshCellSets_.wedges),
            pointOffsets,
            ensightGeometryFile
        );

        writeAllPrims
        (
            "penta6",
            meshCellSets_.nPrisms,
            map(cellShapes, meshCellSets_.prisms),
            pointOffsets,
            ensightGeometryFile
        );

        writeAllPrims
        (
            "pyramid5",
            meshCellSets_.nPyrs,
            map(cellShapes, meshCellSets_.pyrs),
            pointOffsets,
            ensightGeometryFile
        );

        writeAllPrims
        (
            "tetra4",
            meshCellSets_.nTets,
            map(cellShapes, meshCellSets_.tets),
            pointOffsets,
            ensightGeometryFile
        );


        if (meshCellSets_.nPolys)
        {
            if (Pstream::master())
            {
                /*
                ensightGeometryFile
                    << "nfaced" << nl
                    << setw(10) << meshCellSets_.nPolys << nl;
                */
                writePolys
                (
                    meshCellSets_.polys,
                    cellFaces,
                    faces,
                    0,
                    ensightGeometryFile
                );

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    IPstream fromSlave(Pstream::scheduled, slave);
                    labelList polys(fromSlave);
                    cellList cellFaces(fromSlave);
                    faceList faces(fromSlave);

                    writePolys
                    (
                        polys,
                        cellFaces,
                        faces,
                        pointOffsets[slave-1],
                        ensightGeometryFile
                    );
                }
            }
            else
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< meshCellSets_.polys << cellFaces << faces;
            }
        }
    }


    label ensightPatchI = patchPartOffset_;

    forAll(allPatchNames_, patchi)
    {
        const word& patchName = allPatchNames_[patchi];
        const labelList& patchProcessors = allPatchProcs_[patchi];

        if (patchNames_.empty() || patchNames_.found(patchName))
        {
            const nFacePrimitives& nfp = nPatchPrims_.find(patchName)();

            const labelList *trisPtr  = NULL;
            const labelList *quadsPtr = NULL;
            const labelList *polysPtr = NULL;

            const pointField *patchPointsPtr = NULL;
            const faceList *patchFacesPtr = NULL;

            if (mesh_.boundary()[patchi].size())
            {
                const polyPatch& p = mesh_.boundaryMesh()[patchi];

                trisPtr  = &boundaryFaceSets_[patchi].tris;
                quadsPtr = &boundaryFaceSets_[patchi].quads;
                polysPtr = &boundaryFaceSets_[patchi].polys;

                patchPointsPtr = &(p.localPoints());
                patchFacesPtr  = &(p.localFaces());
            }

            const labelList& tris = *trisPtr;
            const labelList& quads = *quadsPtr;
            const labelList& polys = *polysPtr;
            const pointField& patchPoints = *patchPointsPtr;
            const faceList& patchFaces = *patchFacesPtr;

            if (nfp.nTris || nfp.nQuads || nfp.nPolys)
            {
                labelList patchPointOffsets(Pstream::nProcs(), 0);

                if (Pstream::master())
                {
                    ensightGeometryFile
                        << "part" << nl
                        << setw(10) << ensightPatchI++ << nl
                        << patchName << nl
                        << "coordinates" << nl
                        << setw(10) << nfp.nPoints
                        << endl;

                    for (direction d=0; d<vector::nComponents; d++)
                    {
                        if (patchPointsPtr)
                        {
                            writePoints
                            (
                                patchPoints.component(d),
                                ensightGeometryFile
                            );
                        }

                        patchPointOffsets = 0;

                        forAll (patchProcessors, i)
                        {
                            if (patchProcessors[i] != 0)
                            {
                                label slave = patchProcessors[i];
                                IPstream fromSlave(Pstream::scheduled, slave);
                                scalarField patchPointsComponent(fromSlave);

                                writePoints
                                (
                                    patchPointsComponent,
                                    ensightGeometryFile
                                );

                                if (i < Pstream::nProcs()-1)
                                {
                                    patchPointOffsets[i+1] =
                                        patchPointOffsets[i]
                                      + patchPointsComponent.size();
                                }
                            }
                            else
                            {
                                if (i < Pstream::nProcs()-1)
                                {
                                    patchPointOffsets[i+1] =
                                        patchPointOffsets[i]
                                      + patchPoints.size();
                                }
                            }
                        }
                    }
                }
                else if (patchPointsPtr)
                {
                    for (direction d=0; d<vector::nComponents; d++)
                    {
                        OPstream toMaster
                        (
                            Pstream::scheduled,
                            Pstream::masterNo()
                        );
                        toMaster<< patchPoints.component(d);
                    }
                }

                writeAllFacePrims
                (
                    "tria3",
                    tris,
                    nfp.nTris,
                    patchFaces,
                    patchPointOffsets,
                    patchProcessors,
                    ensightGeometryFile
                );

                writeAllFacePrims
                (
                    "quad4",
                    quads,
                    nfp.nQuads,
                    patchFaces,
                    patchPointOffsets,
                    patchProcessors,
                    ensightGeometryFile
                );

                writeAllFacePrims
                (
                    "nsided",
                    polys,
                    nfp.nPolys,
                    patchFaces,
                    patchPointOffsets,
                    patchProcessors,
                    ensightGeometryFile
                );
            }
        }
    }

    if (Pstream::master())
    {
        delete ensightGeometryFilePtr;
    }
}


void Foam::ensightMesh::writeBinary
(
    const fileName& postProcPath,
    const word& prepend,
    const label timeIndex,
    Ostream& ensightCaseFile
) const
{
    //const Time& runTime = mesh.time();
    const pointField& points = mesh_.points();
    const cellList& cellFaces = mesh_.cells();
    const faceList& faces = mesh_.faces();
    const cellShapeList& cellShapes = mesh_.cellShapes();

    word timeFile = prepend;

    if (timeIndex == 0)
    {
        timeFile += "000.";
    }
    else if (mesh_.moving())
    {
        timeFile += itoa(timeIndex) + '.';
    }

    // set the filename of the ensight file
    fileName ensightGeometryFileName = timeFile + "mesh";

    std::ofstream *ensightGeometryFilePtr = NULL;

    if (Pstream::master())
    {
        ensightGeometryFilePtr = new std::ofstream
        (
            (postProcPath/ensightGeometryFileName).c_str(),
            ios_base::out | ios_base::binary | ios_base::trunc
        );
        // Check on file opened?
    }

    std::ofstream& ensightGeometryFile = *ensightGeometryFilePtr;

    if (Pstream::master())
    {
        writeEnsDataBinary("C binary", ensightGeometryFile);
        writeEnsDataBinary("EnSight Geometry File", ensightGeometryFile);
        writeEnsDataBinary("written from OpenFOAM", ensightGeometryFile);
        writeEnsDataBinary("node id assign", ensightGeometryFile);
        writeEnsDataBinary("element id assign", ensightGeometryFile);
    }

    labelList pointOffsets(Pstream::nProcs(), 0);

    if (patchNames_.empty())
    {
        label nPoints = points.size();
        Pstream::gather(nPoints, sumOp<label>());

        if (Pstream::master())
        {
            writeEnsDataBinary("part",ensightGeometryFile);
            writeEnsDataBinary(1,ensightGeometryFile);
            writeEnsDataBinary("internalMesh",ensightGeometryFile);
            writeEnsDataBinary("coordinates",ensightGeometryFile);
            writeEnsDataBinary(nPoints,ensightGeometryFile);

            for (direction d=0; d<vector::nComponents; d++)
            {
                //writePointsBinary(points.component(d), ensightGeometryFile);
                writeEnsDataBinary(points.component(d), ensightGeometryFile);
                pointOffsets[0] = points.size();

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    IPstream fromSlave(Pstream::scheduled, slave);
                    scalarField pointsComponent(fromSlave);
                    //writePointsBinary(pointsComponent, ensightGeometryFile);
                    writeEnsDataBinary(pointsComponent, ensightGeometryFile);
                    pointOffsets[slave] =
                        pointOffsets[slave-1]
                      + pointsComponent.size();
                }
            }
        }
        else
        {
            for (direction d=0; d<vector::nComponents; d++)
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< points.component(d);
            }
        }

        writeAllPrimsBinary
        (
            "hexa8",
            meshCellSets_.nHexesWedges,
            map(cellShapes, meshCellSets_.hexes, meshCellSets_.wedges),
            pointOffsets,
            ensightGeometryFile
        );

        writeAllPrimsBinary
        (
            "penta6",
            meshCellSets_.nPrisms,
            map(cellShapes, meshCellSets_.prisms),
            pointOffsets,
            ensightGeometryFile
        );

        writeAllPrimsBinary
        (
            "pyramid5",
            meshCellSets_.nPyrs,
            map(cellShapes, meshCellSets_.pyrs),
            pointOffsets,
            ensightGeometryFile
        );

        writeAllPrimsBinary
        (
            "tetra4",
            meshCellSets_.nTets,
            map(cellShapes, meshCellSets_.tets),
            pointOffsets,
            ensightGeometryFile
        );

        if (meshCellSets_.nPolys)
        {
            if (Pstream::master())
            {
                /*
                ensightGeometryFile
                    << "nfaced" << nl
                    << setw(10) << meshCellSets_.nPolys << nl;
                */
                writePolysBinary
                (
                    meshCellSets_.polys,
                    cellFaces,
                    faces,
                    0,
                    ensightGeometryFile
                );

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    IPstream fromSlave(Pstream::scheduled, slave);
                    labelList polys(fromSlave);
                    cellList cellFaces(fromSlave);
                    faceList faces(fromSlave);

                    writePolysBinary
                    (
                        polys,
                        cellFaces,
                        faces,
                        pointOffsets[slave-1],
                        ensightGeometryFile
                    );
                }
            }
            else
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< meshCellSets_.polys << cellFaces << faces;
            }
        }

    }

    label ensightPatchI = patchPartOffset_;
    label iCount = 0;

    forAll(allPatchNames_, patchi)
    {
        iCount ++;
        const word& patchName = allPatchNames_[patchi];
        const labelList& patchProcessors = allPatchProcs_[patchi];

        if (patchNames_.empty() || patchNames_.found(patchName))
        {
            const nFacePrimitives& nfp = nPatchPrims_.find(patchName)();

            const labelList *trisPtr = NULL;
            const labelList *quadsPtr = NULL;
            const labelList *polysPtr = NULL;

            const pointField *patchPointsPtr = NULL;
            const faceList *patchFacesPtr = NULL;

            if (mesh_.boundary()[patchi].size())
            {
                const polyPatch& p = mesh_.boundaryMesh()[patchi];

                trisPtr = &boundaryFaceSets_[patchi].tris;
                quadsPtr = &boundaryFaceSets_[patchi].quads;
                polysPtr = &boundaryFaceSets_[patchi].polys;

                patchPointsPtr = &(p.localPoints());
                patchFacesPtr = &(p.localFaces());
            }

            const labelList& tris = *trisPtr;
            const labelList& quads = *quadsPtr;
            const labelList& polys = *polysPtr;
            const pointField& patchPoints = *patchPointsPtr;
            const faceList& patchFaces = *patchFacesPtr;

            if (nfp.nTris || nfp.nQuads || nfp.nPolys)
            {
                labelList patchPointOffsets(Pstream::nProcs(), 0);

                if (Pstream::master())
                {
                    writeEnsDataBinary("part",ensightGeometryFile);
                    writeEnsDataBinary(ensightPatchI++,ensightGeometryFile);
                    //writeEnsDataBinary(patchName.c_str(),ensightGeometryFile);
                    writeEnsDataBinary(patchName.c_str(),ensightGeometryFile);
                    writeEnsDataBinary("coordinates",ensightGeometryFile);
                    writeEnsDataBinary(nfp.nPoints,ensightGeometryFile);

                    for (direction d=0; d<vector::nComponents; d++)
                    {
                        if (patchPointsPtr)
                        {
                            //writePointsBinary
                            writeEnsDataBinary
                            (
                                patchPoints.component(d),
                                ensightGeometryFile
                            );
                        }

                        patchPointOffsets = 0;


                        forAll (patchProcessors, i)
                        {
                            if (patchProcessors[i] != 0)
                            {
                                label slave = patchProcessors[i];
                                IPstream fromSlave(Pstream::scheduled, slave);
                                scalarField patchPointsComponent(fromSlave);

                                //writePointsBinary
                                writeEnsDataBinary
                                (
                                    patchPointsComponent,
                                    ensightGeometryFile
                                );

                                if (i < Pstream::nProcs()-1)
                                {
                                    patchPointOffsets[i+1] =
                                        patchPointOffsets[i]
                                      + patchPointsComponent.size();
                                }
                            }
                            else
                            {
                                if (i < Pstream::nProcs()-1)
                                {
                                    patchPointOffsets[i+1] =
                                        patchPointOffsets[i]
                                      + patchPoints.size();
                                }
                            }
                        }
                    }
                }
                else if (patchPointsPtr)
                {
                    for (direction d=0; d<vector::nComponents; d++)
                    {
                        OPstream toMaster
                        (
                            Pstream::scheduled,
                            Pstream::masterNo()
                        );
                        toMaster<< patchPoints.component(d);
                    }
                }

                writeAllFacePrimsBinary
                (
                    "tria3",
                    tris,
                    nfp.nTris,
                    patchFaces,
                    patchPointOffsets,
                    patchProcessors,
                    ensightGeometryFile
                );

                writeAllFacePrimsBinary
                (
                    "quad4",
                    quads,
                    nfp.nQuads,
                    patchFaces,
                    patchPointOffsets,
                    patchProcessors,
                    ensightGeometryFile
                );

                writeAllFacePrimsBinary
                (
                    "nsided",
                    polys,
                    nfp.nPolys,
                    patchFaces,
                    patchPointOffsets,
                    patchProcessors,
                    ensightGeometryFile
                );
            }
        }
    }


    if (Pstream::master())
    {
        delete ensightGeometryFilePtr;
    }
}

// ************************************************************************* //
