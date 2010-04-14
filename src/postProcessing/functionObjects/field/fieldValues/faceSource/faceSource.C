/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "faceSource.H"
#include "fvMesh.H"
#include "cyclicPolyPatch.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "IOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fieldValues
    {
        defineTypeNameAndDebug(faceSource, 0);
    }

    template<>
    const char* NamedEnum<fieldValues::faceSource::sourceType, 2>::
        names[] = {"faceZone", "patch"};

    const NamedEnum<fieldValues::faceSource::sourceType, 2>
        fieldValues::faceSource::sourceTypeNames_;

    template<>
    const char* NamedEnum<fieldValues::faceSource::operationType, 7>::
        names[] =
        {
            "none", "sum", "areaAverage",
            "areaIntegrate", "weightedAverage", "min", "max"
        };

    const NamedEnum<fieldValues::faceSource::operationType, 7>
        fieldValues::faceSource::operationTypeNames_;

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fieldValues::faceSource::setFaceZoneFaces()
{
    label zoneId = mesh().faceZones().findZoneID(sourceName_);

    if (zoneId < 0)
    {
        FatalErrorIn("faceSource::faceSource::setFaceZoneFaces()")
            << type() << " " << name_ << ": "
            << sourceTypeNames_[source_] << "(" << sourceName_ << "):" << nl
            << "    Unknown face zone name: " << sourceName_
            << ". Valid face zones are: " << mesh().faceZones().names()
            << nl << exit(FatalError);
    }

    const faceZone& fZone = mesh().faceZones()[zoneId];

    faceId_.setSize(fZone.size());
    facePatchId_.setSize(fZone.size());
    faceSign_.setSize(fZone.size());

    label count = 0;
    forAll(fZone, i)
    {
        label faceI = fZone[i];

        label faceId = -1;
        label facePatchId = -1;
        if (mesh().isInternalFace(faceI))
        {
            faceId = faceI;
            facePatchId = -1;
        }
        else
        {
            facePatchId = mesh().boundaryMesh().whichPatch(faceI);
            const polyPatch& pp = mesh().boundaryMesh()[facePatchId];
            if (isA<processorPolyPatch>(pp))
            {
                if (refCast<const processorPolyPatch>(pp).owner())
                {
                    faceId = pp.whichFace(faceI);
                }
                else
                {
                    faceId = -1;
                }
            }
            else if (isA<cyclicPolyPatch>(pp))
            {
                label patchFaceI = faceI - pp.start();
                if (patchFaceI < pp.size()/2)
                {
                    faceId = patchFaceI;
                }
                else
                {
                    faceId = -1;
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceId = faceI - pp.start();
            }
            else
            {
                faceId = -1;
                facePatchId = -1;
            }
        }

        if (faceId >= 0)
        {
            if (fZone.flipMap()[i])
            {
                faceSign_[count] = -1;
            }
            else
            {
                faceSign_[count] = 1;
            }
            faceId_[count] = faceId;
            facePatchId_[count] = facePatchId;
            count++;
        }
    }

    faceId_.setSize(count);
    facePatchId_.setSize(count);
    faceSign_.setSize(count);
    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());

    if (debug)
    {
        Pout<< "Original face zone size = " << fZone.size()
            << ", new size = " << count << endl;
    }
}


void Foam::fieldValues::faceSource::setPatchFaces()
{
    label patchId = mesh().boundaryMesh().findPatchID(sourceName_);

    if (patchId < 0)
    {
        FatalErrorIn("faceSource::constructFaceAddressing()")
            << type() << " " << name_ << ": "
            << sourceTypeNames_[source_] << "(" << sourceName_ << "):" << nl
            << "    Unknown patch name: " << sourceName_
            << ". Valid patch names are: "
            << mesh().boundaryMesh().names() << nl
            << exit(FatalError);
    }

    const polyPatch& pp = mesh().boundaryMesh()[patchId];

    label nFaces = pp.size();
    if (isA<cyclicPolyPatch>(pp))
    {
        nFaces /= 2;
    }
    else if (isA<emptyPolyPatch>(pp))
    {
        nFaces = 0;
    }

    faceId_.setSize(nFaces);
    facePatchId_.setSize(nFaces);
    faceSign_.setSize(nFaces);
    nFaces_ = returnReduce(faceId_.size(), sumOp<label>());

    forAll(faceId_, faceI)
    {
        faceId_[faceI] = faceI;
        facePatchId_[faceI] = patchId;
        faceSign_[faceI] = 1;
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fieldValues::faceSource::initialise(const dictionary& dict)
{
    switch (source_)
    {
        case stFaceZone:
        {
            setFaceZoneFaces();
            break;
        }
        case stPatch:
        {
            setPatchFaces();
            break;
        }
        default:
        {
            FatalErrorIn("faceSource::initialise()")
                << type() << " " << name_ << ": "
                << sourceTypeNames_[source_] << "(" << sourceName_ << "):"
                << nl << "    Unknown source type. Valid source types are:"
                << sourceTypeNames_ << nl << exit(FatalError);
        }
    }

    Info<< type() << " " << name_ << ":" << nl
        << "    total faces  = " << nFaces_
        << nl
        << "    total area   = " << gSum(filterField(mesh().magSf(), false))
        << nl;

    if (operation_ == opWeightedAverage)
    {
        dict.lookup("weightField") >> weightFieldName_;
        if
        (
            obr().foundObject<volScalarField>(weightFieldName_)
         || obr().foundObject<surfaceScalarField>(weightFieldName_)
        )
        {
            Info<< "    weight field = " << weightFieldName_;
        }
        else
        {
            FatalErrorIn("faceSource::initialise()")
                << type() << " " << name_ << ": "
                << sourceTypeNames_[source_] << "(" << sourceName_ << "):"
                << nl << "    Weight field " << weightFieldName_
                << " must be either a " << volScalarField::typeName << " or "
                << surfaceScalarField::typeName << nl << exit(FatalError);
        }
    }

    Info<< nl << endl;
}


void Foam::fieldValues::faceSource::writeFileHeader()
{
    if (outputFilePtr_.valid())
    {
        outputFilePtr_()
            << "# Source : " << sourceTypeNames_[source_] << " "
            << sourceName_ <<  nl << "# Faces  : " << nFaces_ << nl
            << "# Time" << tab << "sum(magSf)";

        forAll(fields_, i)
        {
            outputFilePtr_()
                << tab << operationTypeNames_[operation_]
                << "(" << fields_[i] << ")";
        }

        outputFilePtr_() << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldValues::faceSource::faceSource
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    fieldValue(name, obr, dict, loadFromFiles),
    source_(sourceTypeNames_.read(dict.lookup("source"))),
    operation_(operationTypeNames_.read(dict.lookup("operation"))),
    nFaces_(0),
    faceId_(),
    facePatchId_(),
    faceSign_(),
    weightFieldName_("undefinedWeightedFieldName")
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldValues::faceSource::~faceSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldValues::faceSource::read(const dictionary& dict)
{
    fieldValue::read(dict);

    if (active_)
    {
        initialise(dict);
    }
}


void Foam::fieldValues::faceSource::write()
{
    fieldValue::write();

    if (active_)
    {
        if (Pstream::master())
        {
            outputFilePtr_()
                << obr_.time().value() << tab
                << sum(filterField(mesh().magSf(), false));
        }

        forAll(fields_, i)
        {
            writeValues<scalar>(fields_[i]);
            writeValues<vector>(fields_[i]);
            writeValues<sphericalTensor>(fields_[i]);
            writeValues<symmTensor>(fields_[i]);
            writeValues<tensor>(fields_[i]);
        }

        if (Pstream::master())
        {
            outputFilePtr_()<< endl;
        }

        if (log_)
        {
            Info<< endl;
        }
    }
}


// ************************************************************************* //
