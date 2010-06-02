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

#include "directMappedPatchBase.H"
#include "addToRunTimeSelectionTable.H"
#include "ListListOps.H"
#include "meshSearch.H"
#include "meshTools.H"
#include "OFstream.H"
#include "Random.H"
#include "treeDataFace.H"
#include "indexedOctree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directMappedPatchBase, 0);

    template<>
    const char* NamedEnum<directMappedPatchBase::sampleMode, 3>::names[] =
    {
        "nearestCell",
        "nearestPatchFace",
        "nearestFace"
    };

    const NamedEnum<directMappedPatchBase::sampleMode, 3>
        directMappedPatchBase::sampleModeNames_;


    //- Private class for finding nearest
    //  - point+local index
    //  - sqr(distance)
    //  - processor
    typedef Tuple2<pointIndexHit, Tuple2<scalar, label> > nearInfo;

    class nearestEqOp
    {

    public:

        void operator()(nearInfo& x, const nearInfo& y) const
        {
            if (y.first().hit())
            {
                if (!x.first().hit())
                {
                    x = y;
                }
                else if (y.second().first() < x.second().first())
                {
                    x = y;
                }
            }
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::directMappedPatchBase::collectSamples
(
    pointField& samples,
    labelList& patchFaceProcs,
    labelList& patchFaces,
    pointField& patchFc
) const
{

    // Collect all sample points and the faces they come from.
    List<pointField> globalFc(Pstream::nProcs());
    List<pointField> globalSamples(Pstream::nProcs());
    labelListList globalFaces(Pstream::nProcs());

    globalFc[Pstream::myProcNo()] = patch_.faceCentres();
    globalSamples[Pstream::myProcNo()] = globalFc[Pstream::myProcNo()]+offsets_;
    globalFaces[Pstream::myProcNo()] = identity(patch_.size());

    // Distribute to all processors
    Pstream::gatherList(globalSamples);
    Pstream::scatterList(globalSamples);
    Pstream::gatherList(globalFaces);
    Pstream::scatterList(globalFaces);
    Pstream::gatherList(globalFc);
    Pstream::scatterList(globalFc);

    // Rework into straight list
    samples = ListListOps::combine<pointField>
    (
        globalSamples,
        accessOp<pointField>()
    );
    patchFaces = ListListOps::combine<labelList>
    (
        globalFaces,
        accessOp<labelList>()
    );
    patchFc = ListListOps::combine<pointField>
    (
        globalFc,
        accessOp<pointField>()
    );

    patchFaceProcs.setSize(patchFaces.size());
    labelList nPerProc
    (
        ListListOps::subSizes
        (
            globalFaces,
            accessOp<labelList>()
        )
    );
    label sampleI = 0;
    forAll(nPerProc, procI)
    {
        for (label i = 0; i < nPerProc[procI]; i++)
        {
            patchFaceProcs[sampleI++] = procI;
        }
    }
}


// Find the processor/cell containing the samples. Does not account
// for samples being found in two processors.
void Foam::directMappedPatchBase::findSamples
(
    const pointField& samples,
    labelList& sampleProcs,
    labelList& sampleIndices,
    pointField& sampleLocations
) const
{
    // Lookup the correct region
    const polyMesh& mesh = sampleMesh();

    // All the info for nearest. Construct to miss
    List<nearInfo> nearest(samples.size());

    switch (mode_)
    {
        case NEARESTCELL:
        {
            if (samplePatch_.size() && samplePatch_ != "none")
            {
                FatalErrorIn
                (
                    "directMappedPatchBase::findSamples(const pointField&,"
                    " labelList&, labelList&, pointField&) const"
                )   << "No need to supply a patch name when in "
                    << sampleModeNames_[mode_] << " mode." << exit(FatalError);
            }

            // Octree based search engine
            meshSearch meshSearchEngine(mesh, false);

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                label cellI = meshSearchEngine.findCell(sample);

                if (cellI == -1)
                {
                    nearest[sampleI].second().first() = Foam::sqr(GREAT);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
                else
                {
                    const point& cc = mesh.cellCentres()[cellI];

                    nearest[sampleI].first() = pointIndexHit
                    (
                        true,
                        cc,
                        cellI
                    );
                    nearest[sampleI].second().first() = magSqr(cc-sample);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
            }
            break;
        }

        case NEARESTPATCHFACE:
        {
            Random rndGen(123456);

            const polyPatch& pp = samplePolyPatch();

            if (pp.empty())
            {
                forAll(samples, sampleI)
                {
                    nearest[sampleI].second().first() = Foam::sqr(GREAT);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
            }
            else
            {
                // patch faces
                const labelList patchFaces(identity(pp.size()) + pp.start());

                treeBoundBox patchBb
                (
                    treeBoundBox(pp.points(), pp.meshPoints()).extend
                    (
                        rndGen,
                        1E-4
                    )
                );
                patchBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
                patchBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

                indexedOctree<treeDataFace> boundaryTree
                (
                    treeDataFace    // all information needed to search faces
                    (
                        false,      // do not cache bb
                        mesh,
                        patchFaces  // boundary faces only
                    ),
                    patchBb,        // overall search domain
                    8,              // maxLevel
                    10,             // leafsize
                    3.0             // duplicity
                );

                forAll(samples, sampleI)
                {
                    const point& sample = samples[sampleI];

                    pointIndexHit& nearInfo = nearest[sampleI].first();
                    nearInfo = boundaryTree.findNearest
                    (
                        sample,
                        magSqr(patchBb.span())
                    );

                    if (!nearInfo.hit())
                    {
                        nearest[sampleI].second().first() = Foam::sqr(GREAT);
                        nearest[sampleI].second().second() =
                            Pstream::myProcNo();
                    }
                    else
                    {
                        point fc(pp[nearInfo.index()].centre(pp.points()));
                        nearInfo.setPoint(fc);
                        nearest[sampleI].second().first() = magSqr(fc-sample);
                        nearest[sampleI].second().second() =
                            Pstream::myProcNo();
                    }
                }
            }
            break;
        }

        case NEARESTFACE:
        {
            if (samplePatch_.size() && samplePatch_ != "none")
            {
                FatalErrorIn
                (
                    "directMappedPatchBase::findSamples(const pointField&,"
                    " labelList&, labelList&, pointField&) const"
                )   << "No need to supply a patch name when in "
                    << sampleModeNames_[mode_] << " mode." << exit(FatalError);
            }

            // Octree based search engine
            meshSearch meshSearchEngine(mesh, false);

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                label faceI = meshSearchEngine.findNearestFace(sample);

                if (faceI == -1)
                {
                    nearest[sampleI].second().first() = Foam::sqr(GREAT);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
                else
                {
                    const point& fc = mesh.faceCentres()[faceI];

                    nearest[sampleI].first() = pointIndexHit
                    (
                        true,
                        fc,
                        faceI
                    );
                    nearest[sampleI].second().first() = magSqr(fc-sample);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
            }
            break;
        }

        default:
        {
            FatalErrorIn("directMappedPatchBase::findSamples(..)")
                << "problem." << abort(FatalError);
        }
    }


    // Find nearest.
    Pstream::listCombineGather(nearest, nearestEqOp());
    Pstream::listCombineScatter(nearest);

    if (debug)
    {
        Info<< "directMappedPatchBase::findSamples on mesh " << sampleRegion_
            << " : " << endl;
        forAll(nearest, sampleI)
        {
            label procI = nearest[sampleI].second().second();
            label localI = nearest[sampleI].first().index();

            Info<< "    " << sampleI << " coord:"<< samples[sampleI]
                << " found on processor:" << procI
                << " in local cell/face:" << localI
                << " with cc:" << nearest[sampleI].first().rawPoint() << endl;
        }
    }

    // Check for samples not being found
    forAll(nearest, sampleI)
    {
        if (!nearest[sampleI].first().hit())
        {
            FatalErrorIn
            (
                "directMappedPatchBase::findSamples"
                "(const pointField&, labelList&"
                ", labelList&, pointField&)"
            )   << "Did not find sample " << samples[sampleI]
                << " on any processor of region" << sampleRegion_
                << exit(FatalError);
        }
    }


    // Convert back into proc+local index
    sampleProcs.setSize(samples.size());
    sampleIndices.setSize(samples.size());
    sampleLocations.setSize(samples.size());

    forAll(nearest, sampleI)
    {
        sampleProcs[sampleI] = nearest[sampleI].second().second();
        sampleIndices[sampleI] = nearest[sampleI].first().index();
        sampleLocations[sampleI] = nearest[sampleI].first().hitPoint();
    }
}


void Foam::directMappedPatchBase::calcMapping() const
{
    if (mapPtr_.valid())
    {
        FatalErrorIn("directMappedPatchBase::calcMapping() const")
            << "Mapping already calculated" << exit(FatalError);
    }

    if
    (
        gAverage(mag(offsets_)) <= ROOTVSMALL
     && mode_ == NEARESTPATCHFACE
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
     && samplePatch_ == patch_.name()
    )
    {
        WarningIn("directMappedPatchBase::calcMapping() const")
            << "Invalid offset " << offsets_ << endl
            << "Offset is the vector added to the patch face centres to"
            << " find the patch face supplying the data." << endl
            << "Setting it to " << offsets_
            << " on the same patch, on the same region"
            << " will find the faces themselves which does not make sense"
            << " for anything but testing." << endl
            << "patch_:" << patch_.name() << endl
            << "sampleRegion_:" << sampleRegion_ << endl
            << "mode_:" << sampleModeNames_[mode_] << endl
            << "samplePatch_:" << samplePatch_ << endl
            << "offsets_:" << offsets_ << endl;
    }


    // Get global list of all samples and the processor and face they come from.
    pointField samples;
    labelList patchFaceProcs;
    labelList patchFaces;
    pointField patchFc;
    collectSamples(samples, patchFaceProcs, patchFaces, patchFc);

    // Find processor and cell/face samples are in and actual location.
    labelList sampleProcs;
    labelList sampleIndices;
    pointField sampleLocations;
    findSamples(samples, sampleProcs, sampleIndices, sampleLocations);


    // Now we have all the data we need:
    // - where sample originates from (so destination when mapping):
    //   patchFaces, patchFaceProcs.
    // - cell/face sample is in (so source when mapping)
    //   sampleIndices, sampleProcs.

    //forAll(samples, i)
    //{
    //    Info<< i << " need data in region "
    //        << patch_.boundaryMesh().mesh().name()
    //        << " for proc:" << patchFaceProcs[i]
    //        << " face:" << patchFaces[i]
    //        << " at:" << patchFc[i] << endl
    //        << "Found data in region " << sampleRegion_
    //        << " at proc:" << sampleProcs[i]
    //        << " face:" << sampleIndices[i]
    //        << " at:" << sampleLocations[i]
    //        << nl << endl;
    //}



    if (debug && Pstream::master())
    {
        OFstream str
        (
            patch_.boundaryMesh().mesh().time().path()
          / patch_.name()
          + "_directMapped.obj"
        );
        Pout<< "Dumping mapping as lines from patch faceCentres to"
            << " sampled cell/faceCentres to file " << str.name() << endl;

        label vertI = 0;

        forAll(patchFc, i)
        {
            meshTools::writeOBJ(str, patchFc[i]);
            vertI++;
            meshTools::writeOBJ(str, sampleLocations[i]);
            vertI++;
            str << "l " << vertI-1 << ' ' << vertI << nl;
        }
    }


    //// Check that actual offset vector (sampleLocations - patchFc) is more or
    //// less constant.
    //if (Pstream::master())
    //{
    //    const scalarField magOffset(mag(sampleLocations - patchFc));
    //    const scalar avgOffset(average(magOffset));
    //
    //    forAll(magOffset, sampleI)
    //    {
    //        if
    //        (
    //            mag(magOffset[sampleI]-avgOffset)
    //          > max(SMALL, 0.001*avgOffset)
    //        )
    //        {
    //            WarningIn("directMappedPatchBase::calcMapping() const")
    //                << "The actual cell/face centres picked up using offset "
    //                << offsets_ << " are not" << endl
    //                << "    on a single plane."
    //                << " This might give numerical problems." << endl
    //                << "    At patchface " << patchFc[sampleI]
    //                << " the sampled cell/face " << sampleLocations[sampleI]
    //                << endl
    //                << "    is not on a plane " << avgOffset
    //                << " offset from the patch." << endl
    //                << "    You might want to shift your plane offset."
    //                << " Set the debug flag to get a dump of sampled cells."
    //                << endl;
    //            break;
    //        }
    //    }
    //}


    // Determine schedule.
    mapPtr_.reset(new mapDistribute(sampleProcs, patchFaceProcs));

    // Rework the schedule from indices into samples to cell data to send,
    // face data to receive.

    labelListList& subMap = mapPtr_().subMap();
    labelListList& constructMap = mapPtr_().constructMap();

    forAll(subMap, procI)
    {
        subMap[procI] = UIndirectList<label>
        (
            sampleIndices,
            subMap[procI]
        );
        constructMap[procI] = UIndirectList<label>
        (
            patchFaces,
            constructMap[procI]
        );

        if (debug)
        {
            Pout<< "To proc:" << procI << " sending values of cells/faces:"
                << subMap[procI] << endl;
            Pout<< "From proc:" << procI << " receiving values of patch faces:"
                << constructMap[procI] << endl;
        }
    }

    // Redo constructSize
    mapPtr_().constructSize() = patch_.size();

    if (debug)
    {
        // Check that all elements get a value.
        PackedBoolList used(patch_.size());
        forAll(constructMap, procI)
        {
            const labelList& map = constructMap[procI];

            forAll(map, i)
            {
                label faceI = map[i];

                if (used[faceI] == 0)
                {
                    used[faceI] = 1;
                }
                else
                {
                    FatalErrorIn("directMappedPatchBase::calcMapping() const")
                        << "On patch " << patch_.name()
                        << " patchface " << faceI
                        << " is assigned to more than once."
                        << abort(FatalError);
                }
            }
        }
        forAll(used, faceI)
        {
            if (used[faceI] == 0)
            {
                FatalErrorIn("directMappedPatchBase::calcMapping() const")
                    << "On patch " << patch_.name()
                    << " patchface " << faceI
                    << " is never assigned to."
                    << abort(FatalError);
            }
        }
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::directMappedPatchBase::directMappedPatchBase
(
    const polyPatch& pp
)
:
    patch_(pp),
    sampleRegion_(patch_.boundaryMesh().mesh().name()),
    mode_(NEARESTPATCHFACE),
    samplePatch_("none"),
    uniformOffset_(true),
    offset_(vector::zero),
    offsets_(pp.size(), offset_),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(NULL)
{}


Foam::directMappedPatchBase::directMappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const vectorField& offsets
)
:
    patch_(pp),
    sampleRegion_(sampleRegion),
    mode_(mode),
    samplePatch_(samplePatch),
    uniformOffset_(false),
    offsets_(offsets),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(NULL)
{}


Foam::directMappedPatchBase::directMappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const vector& offset
)
:
    patch_(pp),
    sampleRegion_(sampleRegion),
    mode_(mode),
    samplePatch_(samplePatch),
    uniformOffset_(true),
    offset_(offset),
    offsets_(pp.size(), offset_),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(NULL)
{}


Foam::directMappedPatchBase::directMappedPatchBase
(
    const polyPatch& pp,
    const dictionary& dict
)
:
    patch_(pp),
    sampleRegion_
    (
        dict.lookupOrDefault
        (
            "sampleRegion",
            patch_.boundaryMesh().mesh().name()
        )
    ),
    mode_(sampleModeNames_.read(dict.lookup("sampleMode"))),
    samplePatch_(dict.lookup("samplePatch")),
    uniformOffset_(dict.found("offset")),
    offset_
    (
        uniformOffset_
      ? point(dict.lookup("offset"))
      : vector::zero
    ),
    offsets_
    (
        uniformOffset_
      ? pointField(pp.size(), offset_)
      : dict.lookup("offsets")
    ),
    sameRegion_(sampleRegion_ == patch_.boundaryMesh().mesh().name()),
    mapPtr_(NULL)
{}


Foam::directMappedPatchBase::directMappedPatchBase
(
    const polyPatch& pp,
    const directMappedPatchBase& dmp
)
:
    patch_(pp),
    sampleRegion_(dmp.sampleRegion_),
    mode_(dmp.mode_),
    samplePatch_(dmp.samplePatch_),
    uniformOffset_(dmp.uniformOffset_),
    offset_(dmp.offset_),
    offsets_(dmp.offsets_),
    sameRegion_(dmp.sameRegion_),
    mapPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::directMappedPatchBase::~directMappedPatchBase()
{
    clearOut();
}


void Foam::directMappedPatchBase::clearOut()
{
    mapPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::polyMesh& Foam::directMappedPatchBase::sampleMesh() const
{
    return patch_.boundaryMesh().mesh().time().lookupObject<polyMesh>
    (
        sampleRegion_
    );
}


const Foam::polyPatch& Foam::directMappedPatchBase::samplePolyPatch() const
{
    const polyMesh& nbrMesh = sampleMesh();

    label patchI = nbrMesh.boundaryMesh().findPatchID(samplePatch_);

    if (patchI == -1)
    {
        FatalErrorIn("directMappedPatchBase::samplePolyPatch() ")
            << "Cannot find patch " << samplePatch_
            << " in region " << sampleRegion_ << endl
            << "Valid patches are " << nbrMesh.boundaryMesh().names()
            << exit(FatalError);
    }

    return nbrMesh.boundaryMesh()[patchI];
}


void Foam::directMappedPatchBase::write(Ostream& os) const
{
    os.writeKeyword("sampleMode") << sampleModeNames_[mode_]
        << token::END_STATEMENT << nl;
    os.writeKeyword("sampleRegion") << sampleRegion_
        << token::END_STATEMENT << nl;
    os.writeKeyword("samplePatch") << samplePatch_
        << token::END_STATEMENT << nl;
    if (uniformOffset_)
    {
        os.writeKeyword("offset") << offset_ << token::END_STATEMENT << nl;
    }
    else
    {
        os.writeKeyword("offsets") << offsets_ << token::END_STATEMENT << nl;
    }
}


// ************************************************************************* //
