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

\*----------------------------------------------------------------------------*/

#include "ensightPartCells.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstream.H"
#include "IStringStream.H"
#include "dictionary.H"
#include "cellModeller.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
   defineTypeNameAndDebug(ensightPartCells, 0);
   addToRunTimeSelectionTable(ensightPart, ensightPartCells, istream);
}

Foam::List<Foam::word> Foam::ensightPartCells::elemTypes_
(
    IStringStream
    (
        "(tetra4 pyramid5 penta6 hexa8 nfaced)"
    )()
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightPartCells::classify(const labelList& idList)
{
    // References to cell shape models
    const cellModel& tet   = *(cellModeller::lookup("tet"));
    const cellModel& pyr   = *(cellModeller::lookup("pyr"));
    const cellModel& prism = *(cellModeller::lookup("prism"));
    const cellModel& hex   = *(cellModeller::lookup("hex"));

    const polyMesh& mesh = *meshPtr_;
    const cellShapeList& cellShapes = mesh.cellShapes();

    offset_ = 0;
    size_ = mesh.nCells();

    bool limited = false;
    if (&idList)
    {
        limited = true;
        size_ = idList.size();
    }

    // count the shapes
    label nTet   = 0;
    label nPyr   = 0;
    label nPrism = 0;
    label nHex   = 0;
    label nPoly  = 0;


    // TODO: allow tet-decomposition of polyhedral cells
#if 0
    label nTetDecomp = 0;
    label nPyrDecomp = 0;
#endif

    for (label listI = 0; listI < size_; ++listI)
    {
        label cellId = listI;
        if (limited)
        {
            cellId = idList[listI];
        }

        const cellShape& cellShape = cellShapes[cellId];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == tet)
        {
            nTet++;
        }
        else if (cellModel == pyr)
        {
            nPyr++;
        }
        else if (cellModel == prism)
        {
            nPrism++;
        }
        else if (cellModel == hex)
        {
            nHex++;
        }
        else
        {
            nPoly++;

            // TODO: allow tet-decomposition of polyhedral cells
#if 0
            const cell& cFaces = mesh.cells()[cellI];

            forAll(cFaces, cFaceI)
            {
                const face& f = mesh.faces()[cFaces[cFaceI]];

                label nQuads = 0;
                label nTris = 0;
                f.nTrianglesQuads(mesh.points(), nTris, nQuads);

                nTetDecomp += nTris;
                nPyrDecomp += nQuads;
            }

            nAddCells--;
            nAddPoints++;
#endif
        }
    }


    // we can avoid double looping, but at the cost of allocation
    labelList tetCells(nTet);
    labelList pyramidCells(nPyr);
    labelList prismCells(nPrism);
    labelList hexCells(nHex);
    labelList polyCells(nPoly);

    nTet   = 0,
    nPyr   = 0;
    nPrism = 0;
    nHex   = 0;
    nPoly  = 0;

    // classify the shapes
    for (label listI = 0; listI < size_; ++listI)
    {
        label cellId = listI;
        if (limited)
        {
            cellId = idList[listI];
        }

        const cellShape& cellShape = cellShapes[cellId];
        const cellModel& cellModel = cellShape.model();

        if (cellModel == tet)
        {
            tetCells[nTet++] = cellId;
        }
        else if (cellModel == pyr)
        {
            pyramidCells[nPyr++] = cellId;
        }
        else if (cellModel == prism)
        {
            prismCells[nPrism++] = cellId;
        }
        else if (cellModel == hex)
        {
            hexCells[nHex++] = cellId;
        }
        else
        {
            polyCells[nPoly++] = cellId;

            // TODO: allow tet-decomposition of polyhedral cells
#if 0
            // Mapping from additional point to cell
            addPointCellLabels_[api] = cellId;

            const cell& cFaces = mesh.cells()[cellId];

            forAll(cFaces, cFaceI)
            {
                const face& f = mesh.faces()[cFaces[cFaceI]];

                label nQuads = 0;
                label nTris = 0;
                f.nTrianglesQuads(mesh.points(), nTris, nQuads);

                nTetDecomp += nTris;
                nPyrDecomp += nQuads;
            }

            nAddCells--;
            nAddPoints++;
#endif
        }
    }


    // MUST match with elementTypes
    elemLists_.setSize(elementTypes().size());

    elemLists_[tetra4Elements].transfer( tetCells );
    elemLists_[pyramid5Elements].transfer( pyramidCells );
    elemLists_[penta6Elements].transfer( prismCells );
    elemLists_[hexa8Elements].transfer( hexCells );
    elemLists_[nfacedElements].transfer( polyCells );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightPartCells::ensightPartCells
(
    label partNumber,
    const string& partDescription
)
:
    ensightPart(partNumber, partDescription)
{}


Foam::ensightPartCells::ensightPartCells
(
    label partNumber,
    const polyMesh& pMesh
)
:
    ensightPart(partNumber, "cells", pMesh)
{
    classify();
}


Foam::ensightPartCells::ensightPartCells
(
    label partNumber,
    const polyMesh& pMesh,
    const labelList& idList
)
:
    ensightPart(partNumber, "cells", pMesh)
{
    classify(idList);
}


Foam::ensightPartCells::ensightPartCells
(
    label partNumber,
    const polyMesh& pMesh,
    const cellZone& cZone
)
:
    ensightPart(partNumber, cZone.name(), pMesh)
{
    classify(cZone);
}


Foam::ensightPartCells::ensightPartCells(const ensightPartCells& part)
:
    ensightPart(part)
{}


Foam::ensightPartCells::ensightPartCells(Istream& is)
:
    ensightPart()
{
    reconstruct(is);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightPartCells::~ensightPartCells()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::ensightPart::localPoints Foam::ensightPartCells::calcLocalPoints() const
{
    const polyMesh& mesh = *meshPtr_;

    localPoints ptList(mesh);
    labelList& usedPoints = ptList.list;
    label nPoints = 0;

    forAll(elemLists_, typeI)
    {
        const labelList& idList = elemLists_[typeI];

        // add all points from cells
        forAll(idList, i)
        {
            label id = idList[i] + offset_;
            const labelList& cFaces = mesh.cells()[id];

            forAll(cFaces, cFaceI)
            {
                const face& f = mesh.faces()[cFaces[cFaceI]];

                forAll(f, fp)
                {
                    if (usedPoints[f[fp]] == -1)
                    {
                        usedPoints[f[fp]] = nPoints++;
                    }
                }
            }
        }
    }

    // this is not absolutely necessary, but renumber anyhow
    nPoints = 0;
    forAll(usedPoints, ptI)
    {
        if (usedPoints[ptI] > -1)
        {
            usedPoints[ptI] = nPoints++;
        }
    }

    ptList.nPoints = nPoints;
    return ptList;
}


void Foam::ensightPartCells::writeConnectivity
(
    ensightGeoFile& os,
    const string& key,
    const labelList& idList,
    const labelList& pointMap
) const
{
    os.writeKeyword(key);
    os.write(idList.size());
    os.newline();

    const polyMesh& mesh = *meshPtr_;

    // write polyhedral
    if (word(key) == "nfaced")
    {
        const faceList& meshFaces = mesh.faces();

        // write the number of faces per element
        forAll(idList, i)
        {
            label id = idList[i] + offset_;
            const labelList& cFace = mesh.cells()[id];

            os.write( cFace.size() );
            os.newline();
        }

        // write the number of points per element face
        forAll(idList, i)
        {
            label id = idList[i] + offset_;
            const labelList& cFace = mesh.cells()[id];

            forAll(cFace, faceI)
            {
                const face& cf = meshFaces[cFace[faceI]];

                os.write( cf.size() );
                os.newline();
            }
        }

        // write the points describing each element face
        forAll(idList, i)
        {
            label id = idList[i] + offset_;
            const labelList& cFace = mesh.cells()[id];

            forAll(cFace, faceI)
            {
                const face& cf = meshFaces[cFace[faceI]];

                forAll(cf, ptI)
                {
                    // convert global -> local index
                    // (note: Ensight indices start with 1)
                    os.write( pointMap[cf[ptI]] + 1);
                }
                os.newline();
            }
        }
    }
    else
    {
        // write primitive
        const cellShapeList& cellShapes = mesh.cellShapes();

        forAll(idList, i)
        {
            label id = idList[i] + offset_;
            const cellShape& cellPoints = cellShapes[id];

            // convert global -> local index
            // (note: Ensight indices start with 1)
            forAll(cellPoints, ptI)
            {
                os.write( pointMap[cellPoints[ptI]] + 1 );
            }
            os.newline();
        }
    }
}


// ************************************************************************* //
