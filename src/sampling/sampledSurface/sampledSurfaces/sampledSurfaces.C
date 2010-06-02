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

#include "sampledSurfaces.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IOmanip.H"
#include "ListListOps.H"
#include "mergePoints.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    //- Used to offset faces in Pstream::combineOffset
    template <>
    class offsetOp<face>
    {

    public:

        face operator()
        (
            const face& x,
            const label offset
        ) const
        {
            face result(x.size());

            forAll(x, xI)
            {
                result[xI] = x[xI] + offset;
            }
            return result;
        }
    };


    defineTypeNameAndDebug(sampledSurfaces, 0);
}


bool Foam::sampledSurfaces::verbose_(false);
Foam::scalar Foam::sampledSurfaces::mergeTol_(1e-10);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::sampledSurfaces::classifyFieldTypes()
{
    label nFields = 0;

    scalarFields_.clear();
    vectorFields_.clear();
    sphericalTensorFields_.clear();
    symmTensorFields_.clear();
    tensorFields_.clear();

    forAll(fieldNames_, fieldI)
    {
        const word& fieldName = fieldNames_[fieldI];
        word fieldType = "";

        // check files for a particular time
        if (loadFromFiles_)
        {
            IOobject io
            (
                fieldName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            if (io.headerOk())
            {
                fieldType = io.headerClassName();
            }
            else
            {
                continue;
            }
        }
        else
        {
            // check objectRegistry
            objectRegistry::const_iterator iter = mesh_.find(fieldName);

            if (iter != mesh_.objectRegistry::end())
            {
                fieldType = iter()->type();
            }
            else
            {
                continue;
            }
        }


        if (fieldType == volScalarField::typeName)
        {
            scalarFields_.append(fieldName);
            nFields++;
        }
        else if (fieldType == volVectorField::typeName)
        {
            vectorFields_.append(fieldName);
            nFields++;
        }
        else if (fieldType == volSphericalTensorField::typeName)
        {
            sphericalTensorFields_.append(fieldName);
            nFields++;
        }
        else if (fieldType == volSymmTensorField::typeName)
        {
            symmTensorFields_.append(fieldName);
            nFields++;
        }
        else if (fieldType == volTensorField::typeName)
        {
            tensorFields_.append(fieldName);
            nFields++;
        }

    }

    return nFields;
}


void Foam::sampledSurfaces::writeGeometry() const
{
    // Write to time directory under outputPath_
    // skip surface without faces (eg, a failed cut-plane)

    const fileName outputDir = outputPath_/mesh_.time().timeName();

    forAll(*this, surfI)
    {
        const sampledSurface& s = operator[](surfI);

        if (Pstream::parRun())
        {
            if (Pstream::master() && mergeList_[surfI].faces.size())
            {
                genericFormatter_->write
                (
                    outputDir,
                    s.name(),
                    mergeList_[surfI].points,
                    mergeList_[surfI].faces
                );
            }
        }
        else if (s.faces().size())
        {
            genericFormatter_->write
            (
                outputDir,
                s.name(),
                s.points(),
                s.faces()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::sampledSurfaces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    PtrList<sampledSurface>(),
    name_(name),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    outputPath_(fileName::null),
    fieldNames_(),
    interpolationScheme_(word::null),
    writeFormat_(word::null),
    mergeList_(),
    genericFormatter_(NULL),
    scalarFields_(),
    vectorFields_(),
    sphericalTensorFields_(),
    symmTensorFields_(),
    tensorFields_()
{
    if (Pstream::parRun())
    {
        outputPath_ = mesh_.time().path()/".."/name_;
    }
    else
    {
        outputPath_ = mesh_.time().path()/name_;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledSurfaces::~sampledSurfaces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sampledSurfaces::verbose(const bool verbosity)
{
    verbose_ = verbosity;
}


void Foam::sampledSurfaces::execute()
{
    // Do nothing - only valid on write
}


void Foam::sampledSurfaces::end()
{
    // Do nothing - only valid on write
}


void Foam::sampledSurfaces::write()
{
    if (size())
    {
        // finalize surfaces, merge points etc.
        update();

        const label nFields = classifyFieldTypes();

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

                Pout<< "Creating directory "
                    << outputPath_/mesh_.time().timeName() << nl << endl;

            }

            mkDir(outputPath_/mesh_.time().timeName());
        }

        // write geometry first if required, or when no fields would otherwise
        // be written
        if (nFields == 0 || genericFormatter_->separateFiles())
        {
            writeGeometry();
        }

        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);
    }
}


void Foam::sampledSurfaces::read(const dictionary& dict)
{
    fieldNames_ = wordList(dict.lookup("fields"));

    const label nFields = fieldNames_.size();

    scalarFields_.reset(nFields);
    vectorFields_.reset(nFields);
    sphericalTensorFields_.reset(nFields);
    symmTensorFields_.reset(nFields);
    tensorFields_.reset(nFields);

    interpolationScheme_ = dict.lookupOrDefault<word>
    (
        "interpolationScheme",
        "cell"
    );
    writeFormat_ = dict.lookupOrDefault<word>
    (
        "surfaceFormat",
        "null"
    );


    // define the generic (geometry) writer
    genericFormatter_ = surfaceWriter<bool>::New(writeFormat_);


    PtrList<sampledSurface> newList
    (
        dict.lookup("surfaces"),
        sampledSurface::iNew(mesh_)
    );

    transfer(newList);

    if (Pstream::parRun())
    {
        mergeList_.setSize(size());
    }

    // ensure all surfaces and merge information are expired
    expire();

    if (Pstream::master() && debug)
    {
        Pout<< "sample fields:" << fieldNames_ << nl
            << "sample surfaces:" << nl << "(" << nl;

        forAll(*this, surfI)
        {
            Pout<< "  " << operator[](surfI) << endl;
        }
        Pout<< ")" << endl;
    }
}


void Foam::sampledSurfaces::updateMesh(const mapPolyMesh&)
{
    expire();
}


void Foam::sampledSurfaces::movePoints(const pointField&)
{
    expire();
}


void Foam::sampledSurfaces::readUpdate(const polyMesh::readUpdateState state)
{
    if (state != polyMesh::UNCHANGED)
    {
        expire();
    }
}


bool Foam::sampledSurfaces::needsUpdate() const
{
    forAll(*this, surfI)
    {
        if (operator[](surfI).needsUpdate())
        {
            return true;
        }
    }

    return false;
}


bool Foam::sampledSurfaces::expire()
{
    bool justExpired = false;

    forAll(*this, surfI)
    {
        if (operator[](surfI).expire())
        {
            justExpired = true;
        }

        // clear merge information
        if (Pstream::parRun())
        {
            mergeList_[surfI].clear();
        }
    }

    // reset interpolation
    pointMesh::Delete(mesh_);
    volPointInterpolation::Delete(mesh_);

    // true if any surfaces just expired
    return justExpired;
}


bool Foam::sampledSurfaces::update()
{
    bool updated = false;

    if (!needsUpdate())
    {
        return updated;
    }

    // serial: quick and easy, no merging required
    if (!Pstream::parRun())
    {
        forAll(*this, surfI)
        {
            if (operator[](surfI).update())
            {
                updated = true;
            }
        }

        return updated;
    }

    // dimension as fraction of mesh bounding box
    scalar mergeDim = mergeTol_ * mesh_.globalData().bb().mag();

    if (Pstream::master() && debug)
    {
        Pout<< nl << "Merging all points within "
            << mergeDim << " meter" << endl;
    }

    forAll(*this, surfI)
    {
        sampledSurface& s = operator[](surfI);

        if (s.update())
        {
            updated = true;
        }
        else
        {
            continue;
        }


        // Collect points from all processors
        List<pointField> gatheredPoints(Pstream::nProcs());
        gatheredPoints[Pstream::myProcNo()] = s.points();
        Pstream::gatherList(gatheredPoints);

        if (Pstream::master())
        {
            mergeList_[surfI].points = ListListOps::combine<pointField>
            (
                gatheredPoints,
                accessOp<pointField>()
            );
        }

        // Collect faces from all processors and renumber using sizes of
        // gathered points
        List<faceList> gatheredFaces(Pstream::nProcs());
        gatheredFaces[Pstream::myProcNo()] = s.faces();
        Pstream::gatherList(gatheredFaces);

        if (Pstream::master())
        {
            mergeList_[surfI].faces = static_cast<const faceList&>
            (
                ListListOps::combineOffset<faceList>
                (
                    gatheredFaces,
                    ListListOps::subSizes
                    (
                        gatheredPoints,
                        accessOp<pointField>()
                    ),
                    accessOp<faceList>(),
                    offsetOp<face>()
                )
            );
        }

        pointField newPoints;
        labelList oldToNew;

        bool hasMerged = mergePoints
        (
            mergeList_[surfI].points,
            mergeDim,
            false,                  // verbosity
            oldToNew,
            newPoints
        );

        if (hasMerged)
        {
            // Store point mapping
            mergeList_[surfI].pointsMap.transfer(oldToNew);

            // Copy points
            mergeList_[surfI].points.transfer(newPoints);

            // Relabel faces
            faceList& faces = mergeList_[surfI].faces;

            forAll(faces, faceI)
            {
                inplaceRenumber(mergeList_[surfI].pointsMap, faces[faceI]);
            }

            if (Pstream::master() && debug)
            {
                Pout<< "For surface " << surfI << " merged from "
                    << mergeList_[surfI].pointsMap.size() << " points down to "
                    << mergeList_[surfI].points.size()    << " points" << endl;
            }
        }
    }

    return updated;
}


// ************************************************************************* //
