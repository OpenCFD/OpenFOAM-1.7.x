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

#include "ensightField.H"
#include "fvMesh.H"
#include "volFields.H"
#include "OFstream.H"
#include "IOmanip.H"
#include "itoa.H"
#include "ensightWriteBinary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void writeData(const scalarField& sf, OFstream& ensightFile)
{
    forAll(sf, i)
    {
        if (mag( sf[i] ) >= scalar(floatScalarVSMALL))
        {
            ensightFile << setw(12) << sf[i] << nl;
        }
        else
        {
            ensightFile << setw(12) << scalar(0) << nl;
        }
    }
}


template<class Type>
scalarField map
(
    const Field<Type>& vf,
    const labelList& map,
    const label cmpt
)
{
    scalarField mf(map.size());

    forAll(map, i)
    {
        mf[i] = component(vf[map[i]], cmpt);
    }

    return mf;
}


template<class Type>
scalarField map
(
    const Field<Type>& vf,
    const labelList& map1,
    const labelList& map2,
    const label cmpt
)
{
    scalarField mf(map1.size() + map2.size());

    forAll(map1, i)
    {
        mf[i] = component(vf[map1[i]], cmpt);
    }

    label offset = map1.size();

    forAll(map2, i)
    {
        mf[i + offset] = component(vf[map2[i]], cmpt);
    }

    return mf;
}


template<class Type>
void writeAllData
(
    const char* key,
    const Field<Type>& vf,
    const labelList& prims,
    const label nPrims,
    OFstream& ensightFile
)
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            ensightFile << key << nl;

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                writeData(map(vf, prims, cmpt), ensightFile);

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    IPstream fromSlave(Pstream::scheduled, slave);
                    scalarField data(fromSlave);
                    writeData(data, ensightFile);
                }
            }
        }
        else
        {
            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< map(vf, prims, cmpt);
            }
        }
    }
}


template<class Type>
void writeAllDataBinary
(
    const char* key,
    const Field<Type>& vf,
    const labelList& prims,
    const label nPrims,
    std::ofstream& ensightFile
)
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            writeEnsDataBinary(key,ensightFile);

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                writeEnsDataBinary(map(vf, prims, cmpt), ensightFile);

                for (int slave=1; slave<Pstream::nProcs(); slave++)
                {
                    IPstream fromSlave(Pstream::scheduled, slave);
                    scalarField data(fromSlave);
                    writeEnsDataBinary(data, ensightFile);
                }
            }
        }
        else
        {
            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< map(vf, prims, cmpt);
            }
        }
    }
}


template<class Type>
void writeAllFaceData
(
    const char* key,
    const labelList& prims,
    const label nPrims,
    const Field<Type>& pf,
    const labelList& patchProcessors,
    OFstream& ensightFile
)
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            ensightFile << key << nl;

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                writeData(map(pf, prims, cmpt), ensightFile);

                forAll (patchProcessors, i)
                {
                    if (patchProcessors[i] != 0)
                    {
                        label slave = patchProcessors[i];
                        IPstream fromSlave(Pstream::scheduled, slave);
                        scalarField pf(fromSlave);

                        writeData(pf, ensightFile);
                    }
                }
            }
        }
        else
        {
            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< map(pf, prims, cmpt);
            }
        }
    }
}


template<class Type>
void writeAllFaceDataBinary
(
    const char* key,
    const labelList& prims,
    const label nPrims,
    const Field<Type>& pf,
    const labelList& patchProcessors,
    std::ofstream& ensightFile
)
{
    if (nPrims)
    {
        if (Pstream::master())
        {
            writeEnsDataBinary(key,ensightFile);

            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                writeEnsDataBinary(map(pf, prims, cmpt), ensightFile);

                forAll (patchProcessors, i)
                {
                    if (patchProcessors[i] != 0)
                    {
                        label slave = patchProcessors[i];
                        IPstream fromSlave(Pstream::scheduled, slave);
                        scalarField pf(fromSlave);

                        writeEnsDataBinary(pf, ensightFile);
                    }
                }
            }
        }
        else
        {
            for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
            {
                OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                toMaster<< map(pf, prims, cmpt);
            }
        }
    }
}


template<class Type>
bool writePatchField
(
    const Foam::Field<Type>& pf,
    const Foam::label patchi,
    const Foam::label ensightPatchI,
    const Foam::faceSets& boundaryFaceSet,
    const Foam::ensightMesh::nFacePrimitives& nfp,
    const Foam::labelList& patchProcessors,
    Foam::OFstream& ensightFile
)
{
    if (nfp.nTris || nfp.nQuads || nfp.nPolys)
    {
        if (Pstream::master())
        {
            ensightFile
                << "part" << nl
                << setw(10) << ensightPatchI << nl;
        }

        writeAllFaceData
        (
            "tria3",
            boundaryFaceSet.tris,
            nfp.nTris,
            pf,
            patchProcessors,
            ensightFile
        );

        writeAllFaceData
        (
            "quad4",
            boundaryFaceSet.quads,
            nfp.nQuads,
            pf,
            patchProcessors,
            ensightFile
        );

        writeAllFaceData
        (
            "nsided",
            boundaryFaceSet.polys,
            nfp.nPolys,
            pf,
            patchProcessors,
            ensightFile
        );

        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
bool writePatchFieldBinary
(
    const Foam::Field<Type>& pf,
    const Foam::label patchi,
    const Foam::label ensightPatchI,
    const Foam::faceSets& boundaryFaceSet,
    const Foam::ensightMesh::nFacePrimitives& nfp,
    const Foam::labelList& patchProcessors,
    std::ofstream& ensightFile
)
{
    if (nfp.nTris || nfp.nQuads || nfp.nPolys)
    {
        if (Pstream::master())
        {
            writeEnsDataBinary("part",ensightFile);
            writeEnsDataBinary(ensightPatchI,ensightFile);
        }

        writeAllFaceDataBinary
        (
            "tria3",
            boundaryFaceSet.tris,
            nfp.nTris,
            pf,
            patchProcessors,
            ensightFile
        );

        writeAllFaceDataBinary
        (
            "quad4",
            boundaryFaceSet.quads,
            nfp.nQuads,
            pf,
            patchProcessors,
            ensightFile
        );

        writeAllFaceDataBinary
        (
            "nsided",
            boundaryFaceSet.polys,
            nfp.nPolys,
            pf,
            patchProcessors,
            ensightFile
        );

        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
void writePatchField
(
    const Foam::word& fieldName,
    const Foam::Field<Type>& pf,
    const Foam::word& patchName,
    const Foam::ensightMesh& eMesh,
    const Foam::fileName& postProcPath,
    const Foam::word& prepend,
    const Foam::label timeIndex,
    Foam::Ostream& ensightCaseFile
)
{
    const Time& runTime = eMesh.mesh().time();

    const List<faceSets>& boundaryFaceSets = eMesh.boundaryFaceSets();
    const wordList& allPatchNames = eMesh.allPatchNames();
    const List<labelList>& allPatchProcs = eMesh.allPatchProcs();
    const HashTable<ensightMesh::nFacePrimitives>&
        nPatchPrims = eMesh.nPatchPrims();

    label ensightPatchI = eMesh.patchPartOffset();

    label patchi = -1;

    forAll(allPatchNames, i)
    {
        if (allPatchNames[i] == patchName)
        {
            patchi = i;
            break;
        }
        ensightPatchI++;
    }


    const labelList& patchProcessors = allPatchProcs[patchi];

    word pfName = patchName + '.' + fieldName;

    word timeFile = prepend + itoa(timeIndex);

    OFstream *ensightFilePtr = NULL;
    if (Pstream::master())
    {
        if (timeIndex == 0)
        {
            ensightCaseFile.setf(ios_base::left);

            ensightCaseFile
                << pTraits<Type>::typeName
                << " per element:            1       "
                << setw(15) << pfName
                << (' ' + prepend + "***." + pfName).c_str()
                << nl;
        }

        // set the filename of the ensight file
        fileName ensightFileName(timeFile + "." + pfName);
        ensightFilePtr = new OFstream
        (
            postProcPath/ensightFileName,
            runTime.writeFormat(),
            runTime.writeVersion(),
            runTime.writeCompression()
        );
    }

    OFstream& ensightFile = *ensightFilePtr;

    if (Pstream::master())
    {
        ensightFile << pTraits<Type>::typeName << nl;
    }

    if (patchi >= 0)
    {
        writePatchField
        (
            pf,
            patchi,
            ensightPatchI,
            boundaryFaceSets[patchi],
            nPatchPrims.find(patchName)(),
            patchProcessors,
            ensightFile
        );
    }
    else
    {
        faceSets nullFaceSets;

        writePatchField
        (
            Field<Type>(),
            -1,
            ensightPatchI,
            nullFaceSets,
            nPatchPrims.find(patchName)(),
            patchProcessors,
            ensightFile
        );
    }

    if (Pstream::master())
    {
        delete ensightFilePtr;
    }
}


template<class Type>
void ensightFieldAscii
(
    const Foam::IOobject& fieldObject,
    const Foam::ensightMesh& eMesh,
    const Foam::fileName& postProcPath,
    const Foam::word& prepend,
    const Foam::label timeIndex,
    Foam::Ostream& ensightCaseFile
)
{
    Info<< "Converting field " << fieldObject.name() << endl;

    word timeFile = prepend + itoa(timeIndex);

    const fvMesh& mesh = eMesh.mesh();
    const Time& runTime = mesh.time();

    const cellSets& meshCellSets = eMesh.meshCellSets();
    const List<faceSets>& boundaryFaceSets = eMesh.boundaryFaceSets();
    const wordList& allPatchNames = eMesh.allPatchNames();
    const List<labelList>& allPatchProcs = eMesh.allPatchProcs();
    const wordHashSet& patchNames = eMesh.patchNames();
    const HashTable<ensightMesh::nFacePrimitives>&
        nPatchPrims = eMesh.nPatchPrims();

    const labelList& tets = meshCellSets.tets;
    const labelList& pyrs = meshCellSets.pyrs;
    const labelList& prisms = meshCellSets.prisms;
    const labelList& wedges = meshCellSets.wedges;
    const labelList& hexes = meshCellSets.hexes;
    const labelList& polys = meshCellSets.polys;

    OFstream *ensightFilePtr = NULL;
    if (Pstream::master())
    {
        // set the filename of the ensight file
        fileName ensightFileName(timeFile + "." + fieldObject.name());
        ensightFilePtr = new OFstream
        (
            postProcPath/ensightFileName,
            runTime.writeFormat(),
            runTime.writeVersion(),
            runTime.writeCompression()
        );
    }

    OFstream& ensightFile = *ensightFilePtr;

    GeometricField<Type, fvPatchField, volMesh> vf(fieldObject, mesh);

    if (patchNames.empty())
    {
        if (Pstream::master())
        {
            if (timeIndex == 0)
            {
                ensightCaseFile.setf(ios_base::left);

                ensightCaseFile
                    << pTraits<Type>::typeName
                    << " per element:            1       "
                    << setw(15) << vf.name()
                    << (' ' + prepend + "***." + vf.name()).c_str()
                    << nl;
            }

            ensightFile
                << pTraits<Type>::typeName << nl
                << "part" << nl
                << setw(10) << 1 << nl;

            ensightFile.setf(ios_base::scientific, ios_base::floatfield);
            ensightFile.precision(5);
        }

        if (meshCellSets.nHexesWedges)
        {
            if (Pstream::master())
            {
                ensightFile << "hexa8" << nl;

                for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
                {
                    writeData
                    (
                        map(vf, hexes, wedges, cmpt),
                        ensightFile
                    );

                    for (int slave=1; slave<Pstream::nProcs(); slave++)
                    {
                        IPstream fromSlave(Pstream::scheduled, slave);
                        scalarField data(fromSlave);
                        writeData(data, ensightFile);
                    }
                }
            }
            else
            {
                for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
                {
                    OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                    toMaster<< map(vf, hexes, wedges, cmpt);
                }
            }
        }

        writeAllData("penta6", vf, prisms, meshCellSets.nPrisms, ensightFile);
        writeAllData("pyramid5", vf, pyrs, meshCellSets.nPyrs, ensightFile);
        writeAllData("tetra4", vf, tets, meshCellSets.nTets, ensightFile);
        writeAllData("nfaced", vf, polys, meshCellSets.nPolys, ensightFile);
    }

    label ensightPatchI = eMesh.patchPartOffset();

    forAll(allPatchNames, patchi)
    {
        const word& patchName = allPatchNames[patchi];
        const labelList& patchProcessors = allPatchProcs[patchi];

        if (patchNames.empty() || patchNames.found(patchName))
        {
            if (mesh.boundary()[patchi].size())
            {
                if
                (
                    writePatchField
                    (
                        vf.boundaryField()[patchi],
                        patchi,
                        ensightPatchI,
                        boundaryFaceSets[patchi],
                        nPatchPrims.find(patchName)(),
                        patchProcessors,
                        ensightFile
                    )
                )
                {
                    ensightPatchI++;
                }

            }
            else if (Pstream::master())
            {
                faceSets nullFaceSet;

                if
                (
                    writePatchField
                    (
                        Field<Type>(),
                        -1,
                        ensightPatchI,
                        nullFaceSet,
                        nPatchPrims.find(patchName)(),
                        patchProcessors,
                        ensightFile
                    )
                )
                {
                    ensightPatchI++;
                }
            }
        }
    }

    if (Pstream::master())
    {
        delete ensightFilePtr;
    }
}


template<class Type>
void ensightFieldBinary
(
    const Foam::IOobject& fieldObject,
    const Foam::ensightMesh& eMesh,
    const Foam::fileName& postProcPath,
    const Foam::word& prepend,
    const Foam::label timeIndex,
    Foam::Ostream& ensightCaseFile
)
{
    Info<< "Converting field (binary) " << fieldObject.name() << endl;

    word timeFile = prepend + itoa(timeIndex);

    const fvMesh& mesh = eMesh.mesh();
    //const Time& runTime = mesh.time();

    const cellSets& meshCellSets = eMesh.meshCellSets();
    const List<faceSets>& boundaryFaceSets = eMesh.boundaryFaceSets();
    const wordList& allPatchNames = eMesh.allPatchNames();
    const List<labelList>& allPatchProcs = eMesh.allPatchProcs();
    const wordHashSet& patchNames = eMesh.patchNames();
    const HashTable<ensightMesh::nFacePrimitives>&
        nPatchPrims = eMesh.nPatchPrims();

    const labelList& tets = meshCellSets.tets;
    const labelList& pyrs = meshCellSets.pyrs;
    const labelList& prisms = meshCellSets.prisms;
    const labelList& wedges = meshCellSets.wedges;
    const labelList& hexes = meshCellSets.hexes;
    const labelList& polys = meshCellSets.polys;

    std::ofstream *ensightFilePtr = NULL;
    if (Pstream::master())
    {
        // set the filename of the ensight file
        fileName ensightFileName(timeFile + "." + fieldObject.name());
        ensightFilePtr = new std::ofstream
        (
            (postProcPath/ensightFileName).c_str(),
            ios_base::out | ios_base::binary | ios_base::trunc
        );
        // Check on file opened?
    }

    std::ofstream& ensightFile = *ensightFilePtr;

    GeometricField<Type, fvPatchField, volMesh> vf(fieldObject, mesh);

    if (patchNames.empty())
    {
        if (Pstream::master())
        {
            if (timeIndex == 0)
            {
                ensightCaseFile.setf(ios_base::left);

                ensightCaseFile
                    << pTraits<Type>::typeName
                    << " per element:            1       "
                    << setw(15) << vf.name()
                    << (' ' + prepend + "***." + vf.name()).c_str()
                    << nl;
            }

            writeEnsDataBinary(pTraits<Type>::typeName,ensightFile);
            writeEnsDataBinary("part",ensightFile);
            writeEnsDataBinary(1,ensightFile);
        }

        if (meshCellSets.nHexesWedges)
        {
            if (Pstream::master())
            {
                writeEnsDataBinary("hexa8",ensightFile);

                for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
                {
                    writeEnsDataBinary
                    (
                        map(vf, hexes, wedges, cmpt),
                        ensightFile
                    );

                    for (int slave=1; slave<Pstream::nProcs(); slave++)
                    {
                        IPstream fromSlave(Pstream::scheduled, slave);
                        scalarField data(fromSlave);
                        writeEnsDataBinary(data, ensightFile);
                    }
                }
            }
            else
            {
                for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
                {
                    OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
                    toMaster<< map(vf, hexes, wedges, cmpt);
                }
            }
        }

        writeAllDataBinary
        (
            "penta6",
            vf,
            prisms,
            meshCellSets.nPrisms,
            ensightFile
        );

        writeAllDataBinary
        (
            "pyramid5",
            vf,
            pyrs,
            meshCellSets.nPyrs,
            ensightFile
        );

        writeAllDataBinary
        (
            "tetra4",
            vf,
            tets,
            meshCellSets.nTets,
            ensightFile
        );

        writeAllDataBinary
        (
            "nfaced",
            vf,
            polys,
            meshCellSets.nPolys,
            ensightFile
        );
    }

    label ensightPatchI = eMesh.patchPartOffset();

    forAll(allPatchNames, patchi)
    {
        const word& patchName = allPatchNames[patchi];
        const labelList& patchProcessors = allPatchProcs[patchi];

        if (patchNames.empty() || patchNames.found(patchName))
        {
            if (mesh.boundary()[patchi].size())
            {
                if
                (
                    writePatchFieldBinary
                    (
                        vf.boundaryField()[patchi],
                        patchi,
                        ensightPatchI,
                        boundaryFaceSets[patchi],
                        nPatchPrims.find(patchName)(),
                        patchProcessors,
                        ensightFile
                    )
                )
                {
                    ensightPatchI++;
                }

            }
            else if (Pstream::master())
            {
                faceSets nullFaceSet;

                if
                (
                    writePatchFieldBinary
                    (
                        Field<Type>(),
                        -1,
                        ensightPatchI,
                        nullFaceSet,
                        nPatchPrims.find(patchName)(),
                        patchProcessors,
                        ensightFile
                    )
                )
                {
                    ensightPatchI++;
                }
            }
        }
    }

    if (Pstream::master())
    {
        ensightFile.close();
    }
}


template<class Type>
void ensightField
(
    const Foam::IOobject& fieldObject,
    const Foam::ensightMesh& eMesh,
    const Foam::fileName& postProcPath,
    const Foam::word& prepend,
    const Foam::label timeIndex,
    const bool binary,
    Foam::Ostream& ensightCaseFile
)
{
    if (binary)
    {
        ensightFieldBinary<Type>
        (
            fieldObject,
            eMesh,
            postProcPath,
            prepend,
            timeIndex,
            ensightCaseFile
        );
    }
    else
    {
        ensightFieldAscii<Type>
        (
            fieldObject,
            eMesh,
            postProcPath,
            prepend,
            timeIndex,
            ensightCaseFile
        );
    }
}


// ************************************************************************* //
