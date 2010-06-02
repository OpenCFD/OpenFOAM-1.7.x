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

\*---------------------------------------------------------------------------*/

#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "processorPolyPatch.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::polyBoundaryMesh, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::polyBoundaryMesh::ident(const label len)
{
    labelList elems(len);
    forAll(elems, elemI)
    {
        elems[elemI] = elemI;
    }
    return elems;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& mesh
)
:
    polyPatchList(),
    regIOobject(io),
    mesh_(mesh),
    neighbourEdgesPtr_(NULL)
{
    if (readOpt() == IOobject::MUST_READ)
    {
        polyPatchList& patches = *this;

        // Read polyPatchList
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        patches.setSize(patchEntries.size());

        forAll(patches, patchI)
        {
            patches.set
            (
                patchI,
                polyPatch::New
                (
                    patchEntries[patchI].keyword(),
                    patchEntries[patchI].dict(),
                    patchI,
                    *this
                )
            );
        }

        // Check state of IOstream
        is.check
        (
            "polyBoundaryMesh::polyBoundaryMesh"
            "(const IOobject&, const polyMesh&)"
        );

        close();
    }
}


Foam::polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& pm,
    const label size
)
:
    polyPatchList(size),
    regIOobject(io),
    mesh_(pm),
    neighbourEdgesPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::polyBoundaryMesh::~polyBoundaryMesh()
{
    deleteDemandDrivenData(neighbourEdgesPtr_);
}


void Foam::polyBoundaryMesh::clearGeom()
{
    forAll (*this, patchi)
    {
        operator[](patchi).clearGeom();
    }
}


void Foam::polyBoundaryMesh::clearAddressing()
{
    deleteDemandDrivenData(neighbourEdgesPtr_);

    forAll (*this, patchi)
    {
        operator[](patchi).clearAddressing();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyBoundaryMesh::calcGeometry()
{
    forAll(*this, patchi)
    {
        operator[](patchi).initGeometry();
    }

    forAll(*this, patchi)
    {
        operator[](patchi).calcGeometry();
    }
}


const Foam::List<Foam::labelPairList>&
Foam::polyBoundaryMesh::neighbourEdges() const
{
    if (Pstream::parRun())
    {
        WarningIn("polyBoundaryMesh::neighbourEdges() const")
            << "Neighbour edge addressing not correct across parallel"
            << " boundaries." << endl;
    }

    if (!neighbourEdgesPtr_)
    {
        neighbourEdgesPtr_ = new List<labelPairList>(size());
        List<labelPairList>& neighbourEdges = *neighbourEdgesPtr_;

        // Initialize.
        label nEdgePairs = 0;
        forAll(*this, patchi)
        {
            const polyPatch& pp = operator[](patchi);

            neighbourEdges[patchi].setSize(pp.nEdges() - pp.nInternalEdges());

            forAll(neighbourEdges[patchi], i)
            {
                labelPair& edgeInfo = neighbourEdges[patchi][i];

                edgeInfo[0] = -1;
                edgeInfo[1] = -1;
            }

            nEdgePairs += pp.nEdges() - pp.nInternalEdges();
        }

        // From mesh edge (expressed as a point pair so as not to construct
        // point addressing) to patch + relative edge index.
        HashTable<labelPair, edge, Hash<edge> > pointsToEdge(nEdgePairs);

        forAll(*this, patchi)
        {
            const polyPatch& pp = operator[](patchi);

            const edgeList& edges = pp.edges();

            for
            (
                label edgei = pp.nInternalEdges();
                edgei < edges.size();
                edgei++
            )
            {
                // Edge in patch local points
                const edge& e = edges[edgei];

                // Edge in mesh points.
                edge meshEdge(pp.meshPoints()[e[0]], pp.meshPoints()[e[1]]);

                HashTable<labelPair, edge, Hash<edge> >::iterator fnd =
                    pointsToEdge.find(meshEdge);

                if (fnd == pointsToEdge.end())
                {
                    // First occurrence of mesh edge. Store patch and my
                    // local index.
                    pointsToEdge.insert
                    (
                        meshEdge,
                        labelPair
                        (
                            patchi,
                            edgei - pp.nInternalEdges()
                        )
                    );
                }
                else
                {
                    // Second occurrence. Store.
                    const labelPair& edgeInfo = fnd();

                    neighbourEdges[patchi][edgei - pp.nInternalEdges()] =
                        edgeInfo;

                    neighbourEdges[edgeInfo[0]][edgeInfo[1]]
                         = labelPair(patchi, edgei - pp.nInternalEdges());

                    // Found all two occurrences of this edge so remove from
                    // hash to save space. Note that this will give lots of
                    // problems if the polyBoundaryMesh is multiply connected.
                    pointsToEdge.erase(meshEdge);
                }
            }
        }

        if (pointsToEdge.size())
        {
            FatalErrorIn("polyBoundaryMesh::neighbourEdges() const")
                << "Not all boundary edges of patches match up." << nl
                << "Is the outside of your mesh multiply connected?"
                << abort(FatalError);
        }

        forAll(*this, patchi)
        {
            const polyPatch& pp = operator[](patchi);

            const labelPairList& nbrEdges = neighbourEdges[patchi];

            forAll(nbrEdges, i)
            {
                const labelPair& edgeInfo = nbrEdges[i];

                if (edgeInfo[0] == -1 || edgeInfo[1] == -1)
                {
                    label edgeI = pp.nInternalEdges() + i;
                    const edge& e = pp.edges()[edgeI];

                    FatalErrorIn("polyBoundaryMesh::neighbourEdges() const")
                        << "Not all boundary edges of patches match up." << nl
                        << "Edge " << edgeI << " on patch " << pp.name()
                        << " end points " << pp.localPoints()[e[0]] << ' '
                        << pp.localPoints()[e[1]] << " is not matched to an"
                        << " edge on any other patch." << nl
                        << "Is the outside of your mesh multiply connected?"
                        << abort(FatalError);
                }
            }
        }
    }

    return *neighbourEdgesPtr_;
}


Foam::wordList Foam::polyBoundaryMesh::names() const
{
    const polyPatchList& patches = *this;

    wordList t(patches.size());

    forAll (patches, patchI)
    {
        t[patchI] = patches[patchI].name();
    }

    return t;
}


Foam::wordList Foam::polyBoundaryMesh::types() const
{
    const polyPatchList& patches = *this;

    wordList t(patches.size());

    forAll (patches, patchI)
    {
        t[patchI] = patches[patchI].type();
    }

    return t;
}


Foam::wordList Foam::polyBoundaryMesh::physicalTypes() const
{
    const polyPatchList& patches = *this;

    wordList t(patches.size());

    forAll (patches, patchI)
    {
        t[patchI] = patches[patchI].physicalType();
    }

    return t;
}


Foam::label Foam::polyBoundaryMesh::findPatchID(const word& patchName) const
{
    const polyPatchList& patches = *this;

    forAll (patches, patchI)
    {
        if (patches[patchI].name() == patchName)
        {
            return patchI;
        }
    }

    // Patch not found
    if (debug)
    {
        Pout<< "label polyBoundaryMesh::findPatchID(const word& "
            << "patchName) const"
            << "Patch named " << patchName << " not found.  "
            << "List of available patch names: " << names() << endl;
    }

    // Not found, return -1
    return -1;
}


Foam::label Foam::polyBoundaryMesh::whichPatch(const label faceIndex) const
{
    // Find out which patch the current face belongs to by comparing label
    // with patch start labels.
    // If the face is internal, return -1;
    // if it is off the end of the list, abort
    if (faceIndex >= mesh().nFaces())
    {
        FatalErrorIn
        (
            "polyBoundaryMesh::whichPatch(const label faceIndex) const"
        )   << "given label greater than the number of geometric faces"
            << abort(FatalError);
    }

    if (faceIndex < mesh().nInternalFaces())
    {
        return -1;
    }

    forAll (*this, patchI)
    {
        const polyPatch& bp = operator[](patchI);

        if
        (
            faceIndex >= bp.start()
         && faceIndex < bp.start() + bp.size()
        )
        {
            return patchI;
        }
    }

    // If not in any of above, it is trouble!
    FatalErrorIn
    (
        "label polyBoundaryMesh::whichPatch(const label faceIndex) const"
    )   << "Cannot find face " << faceIndex << " in any of the patches "
        << names() << nl
        << "It seems your patches are not consistent with the mesh :"
        << " internalFaces:" << mesh().nInternalFaces()
        << "  total number of faces:" << mesh().nFaces()
        << abort(FatalError);

    return -1;
}


Foam::labelHashSet Foam::polyBoundaryMesh::patchSet
(
    const wordList& patchNames
) const
{
    wordList allPatchNames = names();
    labelHashSet ps(size());

    forAll(patchNames, i)
    {
        // Treat the given patch names as wild-cards and search the set
        // of all patch names for matches
        labelList patchIDs = findStrings(patchNames[i], allPatchNames);

        if (patchIDs.empty())
        {
            WarningIn("polyBoundaryMesh::patchSet(const wordList&)")
                << "Cannot find any patch names matching " << patchNames[i]
                << endl;
        }

        forAll(patchIDs, j)
        {
            ps.insert(patchIDs[j]);
        }
    }

    return ps;
}


bool Foam::polyBoundaryMesh::checkParallelSync(const bool report) const
{
    if (!Pstream::parRun())
    {
        return false;
    }


    const polyBoundaryMesh& bm = *this;

    bool boundaryError = false;

    // Collect non-proc patches and check proc patches are last.
    wordList names(bm.size());
    wordList types(bm.size());

    label nonProcI = 0;

    forAll (bm, patchI)
    {
        if (!isA<processorPolyPatch>(bm[patchI]))
        {
            if (nonProcI != patchI)
            {
                // There is processor patch inbetween normal patches.
                boundaryError = true;

                if (debug || report)
                {
                    Pout<< " ***Problem with boundary patch " << patchI
                        << " named " << bm[patchI].name()
                        << " of type " <<  bm[patchI].type()
                        << ". The patch seems to be preceeded by processor"
                        << " patches. This is can give problems."
                        << endl;
                }
            }
            else
            {
                names[nonProcI] = bm[patchI].name();
                types[nonProcI] = bm[patchI].type();
                nonProcI++;
            }
        }
    }
    names.setSize(nonProcI);
    types.setSize(nonProcI);

    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = names;
    Pstream::gatherList(allNames);
    Pstream::scatterList(allNames);

    List<wordList> allTypes(Pstream::nProcs());
    allTypes[Pstream::myProcNo()] = types;
    Pstream::gatherList(allTypes);
    Pstream::scatterList(allTypes);

    // Have every processor check but only master print error.

    for (label procI = 1; procI < allNames.size(); procI++)
    {
        if
        (
            (allNames[procI] != allNames[0])
         || (allTypes[procI] != allTypes[0])
        )
        {
            boundaryError = true;

            if (debug || (report && Pstream::master()))
            {
                Info<< " ***Inconsistent patches across processors, "
                       "processor 0 has patch names:" << allNames[0]
                    << " patch types:" << allTypes[0]
                    << " processor " << procI << " has patch names:"
                    << allNames[procI]
                    << " patch types:" << allTypes[procI]
                    << endl;
            }
        }
    }

    return boundaryError;
}


bool Foam::polyBoundaryMesh::checkDefinition(const bool report) const
{
    label nextPatchStart = mesh().nInternalFaces();
    const polyBoundaryMesh& bm = *this;

    bool boundaryError = false;

    forAll (bm, patchI)
    {
        if (bm[patchI].start() != nextPatchStart && !boundaryError)
        {
            boundaryError = true;

            Info<< " ****Problem with boundary patch " << patchI
                << " named " << bm[patchI].name()
                << " of type " <<  bm[patchI].type()
                << ". The patch should start on face no " << nextPatchStart
                << " and the patch specifies " << bm[patchI].start()
                << "." << endl
                << "Possibly consecutive patches have this same problem."
                << " Suppressing future warnings." << endl;
        }

        nextPatchStart += bm[patchI].size();
    }

    reduce(boundaryError, orOp<bool>());

    if (boundaryError)
    {
        if (debug || report)
        {
            Pout << " ***Boundary definition is in error." << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info << "    Boundary definition OK." << endl;
        }

        return false;
    }
}


void Foam::polyBoundaryMesh::movePoints(const pointField& p)
{
    polyPatchList& patches = *this;

    forAll(patches, patchi)
    {
        patches[patchi].initMovePoints(p);
    }

    forAll(patches, patchi)
    {
        patches[patchi].movePoints(p);
    }
}


void Foam::polyBoundaryMesh::updateMesh()
{
    deleteDemandDrivenData(neighbourEdgesPtr_);

    polyPatchList& patches = *this;

    forAll(patches, patchi)
    {
        patches[patchi].initUpdateMesh();
    }

    forAll(patches, patchi)
    {
        patches[patchi].updateMesh();
    }
}


void Foam::polyBoundaryMesh::reorder(const UList<label>& oldToNew)
{
    // Change order of patches
    polyPatchList::reorder(oldToNew);

    // Adapt indices
    polyPatchList& patches = *this;

    forAll(patches, patchi)
    {
        patches[patchi].index() = patchi;
    }

    updateMesh();
}


bool Foam::polyBoundaryMesh::writeData(Ostream& os) const
{
    const polyPatchList& patches = *this;

    os  << patches.size() << nl << token::BEGIN_LIST << incrIndent << nl;

    forAll(patches, patchi)
    {
        os  << indent << patches[patchi].name() << nl
            << indent << token::BEGIN_BLOCK << nl
            << incrIndent << patches[patchi] << decrIndent
            << indent << token::END_BLOCK << endl;
    }

    os  << decrIndent << token::END_LIST;

    // Check state of IOstream
    os.check("polyBoundaryMesh::writeData(Ostream& os) const");

    return os.good();
}


bool Foam::polyBoundaryMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    return regIOobject::writeObject(fmt, ver, IOstream::UNCOMPRESSED);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyBoundaryMesh& pbm)
{
    pbm.writeData(os);
    return os;
}


// ************************************************************************* //
