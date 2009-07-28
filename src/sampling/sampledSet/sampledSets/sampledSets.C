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

#include "sampledSets.H"
#include "dictionary.H"
#include "Time.H"
#include "volFields.H"
#include "ListListOps.H"
#include "SortableList.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledSets, 0);
}

bool Foam::sampledSets::verbose_ = false;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::sampledSets::checkFieldTypes()
{
    wordList fieldTypes(fieldNames_.size());

    // check files for a particular time
    if (loadFromFiles_)
    {
        forAll(fieldNames_, fieldi)
        {
            IOobject io
            (
                fieldNames_[fieldi],
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            if (io.headerOk())
            {
                fieldTypes[fieldi] = io.headerClassName();
            }
            else
            {
                fieldTypes[fieldi] = "(notFound)";
            }
        }
    }
    else
    {
        // check objectRegistry
        forAll(fieldNames_, fieldi)
        {
            objectRegistry::const_iterator iter =
                mesh_.find(fieldNames_[fieldi]);

            if (iter != mesh_.objectRegistry::end())
            {
                fieldTypes[fieldi] = iter()->type();
            }
            else
            {
                fieldTypes[fieldi] = "(notFound)";
            }
        }
    }


    label nFields = 0;

    // classify fieldTypes
    nFields += grep(scalarFields_, fieldTypes);
    nFields += grep(vectorFields_, fieldTypes);
    nFields += grep(sphericalTensorFields_, fieldTypes);
    nFields += grep(symmTensorFields_, fieldTypes);
    nFields += grep(tensorFields_, fieldTypes);

    if (Pstream::master())
    {
        if (debug)
        {
            Pout<< "timeName = " << mesh_.time().timeName() << nl
                << "scalarFields    " << scalarFields_ << nl
                << "vectorFields    " << vectorFields_ << nl
                << "sphTensorFields " << sphericalTensorFields_ << nl
                << "symTensorFields " << symmTensorFields_ <<nl
                << "tensorFields    " << tensorFields_ <<nl;
        }

        if (nFields > 0)
        {
            if (debug)
            {
                Pout<< "Creating directory "
                    << outputPath_/mesh_.time().timeName()
                    << nl << endl;
            }

            mkDir(outputPath_/mesh_.time().timeName());
        }
    }

    return nFields > 0;
}


void Foam::sampledSets::combineSampledSets
(
    PtrList<coordSet>& masterSampledSets,
    labelListList& indexSets
)
{
    // Combine sampleSets from processors. Sort by curveDist. Return
    // ordering in indexSets.
    // Note: only master results are valid

    masterSampledSets_.clear();
    masterSampledSets_.setSize(size());
    indexSets_.setSize(size());

    const PtrList<sampledSet>& sampledSets = *this;

    forAll(sampledSets, seti)
    {
        const sampledSet& samplePts = sampledSets[seti];

        // Collect data from all processors
        List<List<point> > gatheredPts(Pstream::nProcs());
        gatheredPts[Pstream::myProcNo()] = samplePts;
        Pstream::gatherList(gatheredPts);

        List<labelList> gatheredSegments(Pstream::nProcs());
        gatheredSegments[Pstream::myProcNo()] = samplePts.segments();
        Pstream::gatherList(gatheredSegments);

        List<scalarList> gatheredDist(Pstream::nProcs());
        gatheredDist[Pstream::myProcNo()] = samplePts.curveDist();
        Pstream::gatherList(gatheredDist);


        // Combine processor lists into one big list.
        List<point> allPts
        (
            ListListOps::combine<List<point> >
            (
                gatheredPts, accessOp<List<point> >()
            )
        );
        labelList allSegments
        (
            ListListOps::combine<labelList>
            (
                gatheredSegments, accessOp<labelList>()
            )
        );
        scalarList allCurveDist
        (
            ListListOps::combine<scalarList>
            (
                gatheredDist, accessOp<scalarList>()
            )
        );

        // Sort curveDist and use to fill masterSamplePts
        SortableList<scalar> sortedDist(allCurveDist);
        indexSets[seti] = sortedDist.indices();

        // Get reference point (note: only master has all points)
        point refPt;

        if (allPts.size())
        {
            refPt = samplePts.getRefPoint(allPts);
        }
        else
        {
            refPt = vector::zero;
        }


        masterSampledSets.set
        (
            seti,
            new coordSet
            (
                samplePts.name(),
                samplePts.axis(),
                UIndirectList<point>(allPts, indexSets[seti]),
                refPt
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSets::sampledSets
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    PtrList<sampledSet>(),
    name_(name),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    outputPath_(fileName::null),
    searchEngine_(mesh_, true),
    fieldNames_(),
    interpolationScheme_(word::null),
    writeFormat_(word::null)
{
    if (Pstream::parRun())
    {
        outputPath_ = mesh_.time().path()/".."/name_;
    }
    else
    {
        outputPath_ = mesh_.time().path()/name_;
    }
    if (mesh_.name() != fvMesh::defaultRegion)
    {
        outputPath_ = outputPath_/mesh_.name();
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSets::~sampledSets()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sampledSets::verbose(const bool verbosity)
{
    verbose_ = verbosity;
}


void Foam::sampledSets::execute()
{
    // Do nothing - only valid on write
}


void Foam::sampledSets::end()
{
    // Do nothing - only valid on write
}


void Foam::sampledSets::write()
{
    if (size() && checkFieldTypes())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);
    }
}


void Foam::sampledSets::read(const dictionary& dict)
{
    dict_ = dict;

    fieldNames_ = wordList(dict_.lookup("fields"));

    interpolationScheme_ = "cell";
    dict_.readIfPresent("interpolationScheme", interpolationScheme_);

    writeFormat_ = "null";
    dict_.readIfPresent("setFormat", writeFormat_);

    scalarFields_.clear();
    vectorFields_.clear();
    sphericalTensorFields_.clear();
    symmTensorFields_.clear();
    tensorFields_.clear();

    PtrList<sampledSet> newList
    (
        dict_.lookup("sets"),
        sampledSet::iNew(mesh_, searchEngine_)
    );
    transfer(newList);
    combineSampledSets(masterSampledSets_, indexSets_);

    if (Pstream::master() && debug)
    {
        Pout<< "sample fields:" << fieldNames_ << nl
            << "sample sets:" << nl << "(" << nl;

        forAll(*this, si)
        {
            Pout << "  " << operator[](si) << endl;
        }
        Pout << ")" << endl;
    }
}


void Foam::sampledSets::correct()
{
    // reset interpolation
    pointMesh::Delete(mesh_);
    volPointInterpolation::Delete(mesh_);

    searchEngine_.correct();

    PtrList<sampledSet> newList
    (
        dict_.lookup("sets"),
        sampledSet::iNew(mesh_, searchEngine_)
    );
    transfer(newList);
    combineSampledSets(masterSampledSets_, indexSets_);
}


void Foam::sampledSets::updateMesh(const mapPolyMesh&)
{
    correct();
}


void Foam::sampledSets::movePoints(const pointField&)
{
    correct();
}


void Foam::sampledSets::readUpdate(const polyMesh::readUpdateState state)
{
    if (state != polyMesh::UNCHANGED)
    {
        correct();
    }
}


// ************************************************************************* //
