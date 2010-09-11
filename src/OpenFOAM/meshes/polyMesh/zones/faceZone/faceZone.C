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
    A subset of mesh faces.

\*---------------------------------------------------------------------------*/

#include "faceZone.H"
#include "addToRunTimeSelectionTable.H"
#include "faceZoneMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"
#include "mapPolyMesh.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZone, 0);
    defineRunTimeSelectionTable(faceZone, dictionary);
    addToRunTimeSelectionTable(faceZone, faceZone, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceZone::calcFaceZonePatch() const
{
    if (debug)
    {
        Info<< "void faceZone::calcFaceZonePatch() const : "
            << "Calculating primitive patch"
            << endl;
    }

    if (patchPtr_)
    {
        FatalErrorIn
        (
            "void faceZone::calcFaceZonePatch() const"
        )   << "primitive face zone patch already calculated"
            << abort(FatalError);
    }

    patchPtr_ =
        new primitiveFacePatch
        (
            faceList(size()),
            zoneMesh().mesh().points()
        );

    primitiveFacePatch& patch = *patchPtr_;

    const faceList& f = zoneMesh().mesh().faces();

    const labelList& addr = *this;
    const boolList& flip = flipMap();

    forAll (addr, faceI)
    {
        if (flip[faceI])
        {
            patch[faceI] = f[addr[faceI]].reverseFace();
        }
        else
        {
            patch[faceI] = f[addr[faceI]];
        }
    }

    if (debug)
    {
        Info<< "void faceZone::calcFaceZonePatch() const : "
            << "Finished calculating primitive patch"
            << endl;
    }
}


const Foam::Map<Foam::label>& Foam::faceZone::faceLookupMap() const
{
    if (!faceLookupMapPtr_)
    {
        calcFaceLookupMap();
    }

    return *faceLookupMapPtr_;
}


void Foam::faceZone::calcFaceLookupMap() const
{
    if (debug)
    {
        Info<< "void faceZone::calcFaceLookupMap() const : "
            << "Calculating face lookup map"
            << endl;
    }

    if (faceLookupMapPtr_)
    {
        FatalErrorIn
        (
            "void faceZone::calcFaceLookupMap() const"
        )   << "face lookup map already calculated"
            << abort(FatalError);
    }

    const labelList& addr = *this;

    faceLookupMapPtr_ = new Map<label>(2*addr.size());
    Map<label>& flm = *faceLookupMapPtr_;

    forAll (addr, faceI)
    {
        flm.insert(addr[faceI], faceI);
    }

    if (debug)
    {
        Info<< "void faceZone::calcFaceLookupMap() const : "
            << "Finished calculating face lookup map"
            << endl;
    }
}


void Foam::faceZone::calcCellLayers() const
{
    if (debug)
    {
        Info<< "void Foam::faceZone::calcCellLayers() const : "
            << "calculating master cells"
            << endl;
    }

    // It is an error to attempt to recalculate edgeCells
    // if the pointer is already set
    if (masterCellsPtr_ || slaveCellsPtr_)
    {
        FatalErrorIn("void faceZone::calcCellLayers() const")
            << "cell layers already calculated"
            << abort(FatalError);
    }
    else
    {
        // Go through all the faces in the master zone.  Choose the
        // master or slave cell based on the face flip

        const labelList& own = zoneMesh().mesh().faceOwner();
        const labelList& nei = zoneMesh().mesh().faceNeighbour();

        const labelList& mf = *this;

        const boolList& faceFlip = flipMap();

        masterCellsPtr_ = new labelList(mf.size());
        labelList& mc = *masterCellsPtr_;

        slaveCellsPtr_ = new labelList(mf.size());
        labelList& sc = *slaveCellsPtr_;

        forAll (mf, faceI)
        {
            label ownCellI = own[mf[faceI]];
            label neiCellI =
            (
                zoneMesh().mesh().isInternalFace(mf[faceI])
              ? nei[mf[faceI]]
              : -1
            );

            if (!faceFlip[faceI])
            {
                // Face is oriented correctly, no flip needed
                mc[faceI] = neiCellI;
                sc[faceI] = ownCellI;
            }
            else
            {
                mc[faceI] = ownCellI;
                sc[faceI] = neiCellI;
            }
        }
        //Info << "masterCells: " << mc << endl;
        //Info << "slaveCells: " << sc << endl;
    }
}


void Foam::faceZone::checkAddressing() const
{
    if (size() != flipMap_.size())
    {
        FatalErrorIn("void Foam::faceZone::checkAddressing() const")
            << "Different sizes of the addressing and flipMap arrays.  "
            << "Size of addressing: " << size()
            << " size of flip map: " << flipMap_.size()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::faceZone::faceZone
(
    const word& name,
    const labelList& addr,
    const boolList& fm,
    const label index,
    const faceZoneMesh& zm
)
:
    labelList(addr),
    name_(name),
    flipMap_(fm),
    index_(index),
    zoneMesh_(zm),
    patchPtr_(NULL),
    masterCellsPtr_(NULL),
    slaveCellsPtr_(NULL),
    mePtr_(NULL),
    faceLookupMapPtr_(NULL)
{
    checkAddressing();
}


Foam::faceZone::faceZone
(
    const word& name,
    const Xfer<labelList>& addr,
    const Xfer<boolList>& fm,
    const label index,
    const faceZoneMesh& zm
)
:
    labelList(addr),
    name_(name),
    flipMap_(fm),
    index_(index),
    zoneMesh_(zm),
    patchPtr_(NULL),
    masterCellsPtr_(NULL),
    slaveCellsPtr_(NULL),
    mePtr_(NULL),
    faceLookupMapPtr_(NULL)
{
    checkAddressing();
}


// Construct from dictionary
Foam::faceZone::faceZone
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faceZoneMesh& zm
)
:
    labelList(dict.lookup("faceLabels")),
    name_(name),
    flipMap_(dict.lookup("flipMap")),
    index_(index),
    zoneMesh_(zm),
    patchPtr_(NULL),
    masterCellsPtr_(NULL),
    slaveCellsPtr_(NULL),
    mePtr_(NULL),
    faceLookupMapPtr_(NULL)
{
    checkAddressing();
}


// Construct given the original zone and resetting the
// face list and zone mesh information
Foam::faceZone::faceZone
(
    const faceZone& fz,
    const labelList& addr,
    const boolList& fm,
    const label index,
    const faceZoneMesh& zm
)
:
    labelList(addr),
    name_(fz.name()),
    flipMap_(fm),
    index_(index),
    zoneMesh_(zm),
    patchPtr_(NULL),
    masterCellsPtr_(NULL),
    slaveCellsPtr_(NULL),
    mePtr_(NULL),
    faceLookupMapPtr_(NULL)
{
    checkAddressing();
}


Foam::faceZone::faceZone
(
    const faceZone& fz,
    const Xfer<labelList>& addr,
    const Xfer<boolList>& fm,
    const label index,
    const faceZoneMesh& zm
)
:
    labelList(addr),
    name_(fz.name()),
    flipMap_(fm),
    index_(index),
    zoneMesh_(zm),
    patchPtr_(NULL),
    masterCellsPtr_(NULL),
    slaveCellsPtr_(NULL),
    mePtr_(NULL),
    faceLookupMapPtr_(NULL)
{
    checkAddressing();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceZone::~faceZone()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::faceZone::whichFace(const label globalFaceID) const
{
    const Map<label>& flm = faceLookupMap();

    Map<label>::const_iterator flmIter = flm.find(globalFaceID);

    if (flmIter == flm.end())
    {
        return -1;
    }
    else
    {
        return flmIter();
    }
}


const Foam::faceZoneMesh& Foam::faceZone::zoneMesh() const
{
    return zoneMesh_;
}


const Foam::primitiveFacePatch& Foam::faceZone::operator()() const
{
    if (!patchPtr_)
    {
        calcFaceZonePatch();
    }

    return *patchPtr_;
}


const Foam::labelList& Foam::faceZone::masterCells() const
{
    if (!masterCellsPtr_)
    {
        calcCellLayers();
    }

    return *masterCellsPtr_;
}


const Foam::labelList& Foam::faceZone::slaveCells() const
{
    if (!slaveCellsPtr_)
    {
        calcCellLayers();
    }

    return *slaveCellsPtr_;
}


const Foam::labelList& Foam::faceZone::meshEdges() const
{
    if (!mePtr_)
    {
        //labelList faceCells(size());
        //
        //const labelList& own = zoneMesh().mesh().faceOwner();
        //
        //const labelList& faceLabels = *this;
        //
        //forAll (faceCells, faceI)
        //{
        //    faceCells[faceI] = own[faceLabels[faceI]];
        //}
        //
        //mePtr_ =
        //    new labelList
        //    (
        //        operator()().meshEdges
        //        (
        //            zoneMesh().mesh().edges(),
        //            zoneMesh().mesh().cellEdges(),
        //            faceCells
        //        )
        //    );

        mePtr_ =
            new labelList
            (
                operator()().meshEdges
                (
                    zoneMesh().mesh().edges(),
                    zoneMesh().mesh().pointEdges()
                )
            );
    }

    return *mePtr_;
}


void Foam::faceZone::clearAddressing()
{
    deleteDemandDrivenData(patchPtr_);

    deleteDemandDrivenData(masterCellsPtr_);
    deleteDemandDrivenData(slaveCellsPtr_);

    deleteDemandDrivenData(mePtr_);
    deleteDemandDrivenData(faceLookupMapPtr_);
}


void Foam::faceZone::resetAddressing
(
    const labelList& addr,
    const boolList& flipMap
)
{
    clearAddressing();
    labelList::operator=(addr);
    flipMap_ = flipMap;
}


void Foam::faceZone::updateMesh(const mapPolyMesh& mpm)
{
    clearAddressing();

    labelList newAddressing(size());
    boolList newFlipMap(flipMap_.size());
    label nFaces = 0;

    const labelList& faceMap = mpm.reverseFaceMap();

    forAll(*this, i)
    {
        label faceI = operator[](i);

        if (faceMap[faceI] >= 0)
        {
            newAddressing[nFaces] = faceMap[faceI];
            newFlipMap[nFaces] = flipMap_[i];       // Keep flip map.
            nFaces++;
        }
    }

    newAddressing.setSize(nFaces);
    newFlipMap.setSize(nFaces);

    transfer(newAddressing);
    flipMap_.transfer(newFlipMap);
}


bool Foam::faceZone::checkDefinition(const bool report) const
{
    const labelList& addr = *this;

    bool boundaryError = false;

    forAll(addr, i)
    {
        if (addr[i] < 0 || addr[i] >= zoneMesh().mesh().faces().size())
        {
            boundaryError = true;

            if (report)
            {
                SeriousErrorIn
                (
                    "bool faceZone::checkDefinition("
                    "const bool report) const"
                )   << "Zone " << name()
                    << " contains invalid face label " << addr[i] << nl
                    << "Valid face labels are 0.."
                    << zoneMesh().mesh().faces().size()-1 << endl;
            }
        }
    }
    return boundaryError;
}


bool Foam::faceZone::checkParallelSync(const bool report) const
{
    const polyMesh& mesh = zoneMesh().mesh();
    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    bool boundaryError = false;


    // Check that zone faces are synced
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        boolList neiZoneFace(mesh.nFaces()-mesh.nInternalFaces(), false);
        boolList neiZoneFlip(mesh.nFaces()-mesh.nInternalFaces(), false);
        forAll(*this, i)
        {
            label faceI = operator[](i);

            if (!mesh.isInternalFace(faceI))
            {
                neiZoneFace[faceI-mesh.nInternalFaces()] = true;
                neiZoneFlip[faceI-mesh.nInternalFaces()] = flipMap()[i];
            }
        }
        boolList myZoneFace(neiZoneFace);
        syncTools::swapBoundaryFaceList(mesh, neiZoneFace, false);
        boolList myZoneFlip(neiZoneFlip);
        syncTools::swapBoundaryFaceList(mesh, neiZoneFlip, false);

        forAll(*this, i)
        {
            label faceI = operator[](i);

            label patchI = bm.whichPatch(faceI);

            if (patchI != -1 && bm[patchI].coupled())
            {
                label bFaceI = faceI-mesh.nInternalFaces();

                // Check face in zone on both sides
                if (myZoneFace[bFaceI] != neiZoneFace[bFaceI])
                {
                    boundaryError = true;

                    if (report)
                    {
                        Pout<< " ***Problem with faceZone " << index()
                            << " named " << name()
                            << ". Face " << faceI
                            << " on coupled patch "
                            << bm[patchI].name()
                            << " is not consistent with its coupled neighbour."
                            << endl;
                    }
                }
                else if (myZoneFlip[bFaceI] == neiZoneFlip[bFaceI])
                {
                    // Flip state should be opposite.
                    boundaryError = true;

                    if (report)
                    {
                        Pout<< " ***Problem with faceZone " << index()
                            << " named " << name()
                            << ". Face " << faceI
                            << " on coupled patch "
                            << bm[patchI].name()
                            << " does not have consistent flipMap"
                            << " across coupled faces."
                            << endl;
                    }
                }
            }
        }
    }

    return returnReduce(boundaryError, orOp<bool>());
}


void Foam::faceZone::movePoints(const pointField& p)
{
    if (patchPtr_)
    {
        patchPtr_->movePoints(p);
    }
}

void Foam::faceZone::write(Ostream& os) const
{
    os  << nl << name()
        << nl << static_cast<const labelList&>(*this)
        << nl << flipMap();
}


void Foam::faceZone::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type() << token::END_STATEMENT << nl;

    writeEntry("faceLabels", os);
    flipMap().writeEntry("flipMap", os);

    os  << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const faceZone& p)
{
    p.write(os);
    os.check("Ostream& operator<<(Ostream& f, const faceZone& p");
    return os;
}


// ************************************************************************* //
