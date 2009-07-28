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

#include "ensightPartFaces.H"
#include "addToRunTimeSelectionTable.H"
#include "IOstreams.H"
#include "IStringStream.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
   defineTypeNameAndDebug(ensightPartFaces, 0);
   addToRunTimeSelectionTable(ensightPart, ensightPartFaces, istream);
}


Foam::List<Foam::word> Foam::ensightPartFaces::elemTypes_
(
    IStringStream
    (
        "(tria3 quad4 nsided)"
    )()
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightPartFaces::ensightPartFaces
(
    label partNumber,
    const string& partDescription
)
:
    ensightPart(partNumber, partDescription)
{
    isCellData_ = false;
}


Foam::ensightPartFaces::ensightPartFaces
(
    label partNumber,
    const polyMesh& pMesh,
    const polyPatch& pPatch
)
:
    ensightPart(partNumber, pPatch.name(), pMesh)
{
    isCellData_ = false;
    offset_ = pPatch.start();
    size_ = pPatch.size();

    // count the shapes
    label nTri  = 0;
    label nQuad = 0;
    label nPoly = 0;

    forAll (pPatch, patchfaceI)
    {
        const face& f = pMesh.faces()[patchfaceI + offset_];

        if (f.size() == 3)
        {
            nTri++;
        }
        else if (f.size() == 4)
        {
            nQuad++;
        }
        else
        {
            nPoly++;
        }
    }

    // we can avoid double looping, but at the cost of allocation

    labelList triCells(nTri);
    labelList quadCells(nQuad);
    labelList polygonCells(nPoly);

    nTri  = 0;
    nQuad = 0;
    nPoly = 0;

    // classify the shapes
    forAll(pPatch, patchfaceI)
    {
        const face& f = pMesh.faces()[patchfaceI + offset_];

        if (f.size() == 3)
        {
            triCells[nTri++] = patchfaceI;
        }
        else if (f.size() == 4)
        {
            quadCells[nQuad++] = patchfaceI;
        }
        else
        {
            polygonCells[nPoly++] = patchfaceI;
        }
    }


    // MUST match with elementTypes
    elemLists_.setSize(elementTypes().size());

    elemLists_[tria3Elements].transfer( triCells );
    elemLists_[quad4Elements].transfer( quadCells );
    elemLists_[nsidedElements].transfer( polygonCells );
}


Foam::ensightPartFaces::ensightPartFaces(const ensightPartFaces &part)
:
    ensightPart(part)
{}


Foam::ensightPartFaces::ensightPartFaces(Istream& is)
:
    ensightPart()
{
    isCellData_ = false;
    reconstruct(is);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightPartFaces::~ensightPartFaces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::ensightPart::localPoints Foam::ensightPartFaces::calcLocalPoints() const
{
    const polyMesh& mesh = *meshPtr_;

    localPoints ptList(mesh);
    labelList& usedPoints = ptList.list;
    label nPoints = 0;

    forAll(elemLists_, typeI)
    {
        const labelList& idList = elemLists_[typeI];

        // add all points from faces
        forAll(idList, i)
        {
            label id = idList[i] + offset_;
            const face& f = mesh.faces()[id];

            forAll(f, fp)
            {
                if (usedPoints[f[fp]] == -1)
                {
                    usedPoints[f[fp]] = nPoints++;
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


void Foam::ensightPartFaces::writeConnectivity
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

    const faceList& meshFaces = meshPtr_->faces();

    // write (polygon) face sizes
    if (word(key) == "nsided")
    {
        // write the number of points per face
        forAll(idList, i)
        {
            label id = idList[i] + offset_;
            const face& f = meshFaces[id];

            os.write( f.size() );
            os.newline();
        }
    }

    // write the points describing the face
    forAll(idList, i)
    {
        label id = idList[i] + offset_;
        const face& f = meshFaces[id];

        // convert global -> local index
        // (note: Ensight indices start with 1)
        forAll(f, fp)
        {
            os.write( pointMap[f[fp]] + 1 );
        }
        os.newline();
    }
}


// ************************************************************************* //
