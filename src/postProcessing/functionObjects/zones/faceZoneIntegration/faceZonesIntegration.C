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

#include "faceZonesIntegration.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IOmanip.H"
#include "ListListOps.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZonesIntegration, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZonesIntegration::faceZonesIntegration
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    log_(false),
    faceZonesSet_(),
    fItems_(),
    faceZonesIntegrationFilePtr_(NULL)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "Foam::faceZonesIntegration::faceZonesIntegration"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating."
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceZonesIntegration::~faceZonesIntegration()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceZonesIntegration::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", false);

        dict.lookup("fields") >> fItems_;

        dict.lookup("faceZones") >> faceZonesSet_;
    }
}


void Foam::faceZonesIntegration::makeFile()
{
    // Create the face Zone file if not already created
    if (faceZonesIntegrationFilePtr_.empty())
    {
        if (debug)
        {
            Info<< "Creating faceZonesIntegration file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName faceZonesIntegrationDir;
            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                faceZonesIntegrationDir =
                    obr_.time().path()/".."/name_/obr_.time().timeName();
            }
            else
            {
                faceZonesIntegrationDir =
                    obr_.time().path()/name_/obr_.time().timeName();
            }

            // Create directory if does not exist.
            mkDir(faceZonesIntegrationDir);

            // Open new file at start up
            faceZonesIntegrationFilePtr_.resize(fItems_.size());

            forAll(fItems_, Ifields)
            {
                const word& fieldName = fItems_[Ifields];

                OFstream* sPtr = new OFstream
                    (
                        faceZonesIntegrationDir/fieldName
                    );

                faceZonesIntegrationFilePtr_.insert(fieldName, sPtr);
            }

            // Add headers to output data
            writeFileHeader();
        }
    }
}


void Foam::faceZonesIntegration::writeFileHeader()
{
    forAllIter(HashPtrTable<OFstream>, faceZonesIntegrationFilePtr_, iter)
    {
        unsigned int w = IOstream::defaultPrecision() + 7;

        OFstream& os = *faceZonesIntegrationFilePtr_[iter.key()];

        os  << "#Time " << setw(w);

        forAll (faceZonesSet_, zoneI)
        {
            const word name = faceZonesSet_[zoneI];
            os  << name << setw(w);
        }

        os  << nl << endl;
    }
}


void Foam::faceZonesIntegration::execute()
{
    // Do nothing - only valid on write
}


void Foam::faceZonesIntegration::end()
{
    // Do nothing - only valid on write
}


void Foam::faceZonesIntegration::write()
{
    if (active_)
    {
        makeFile();

        forAll(fItems_, fieldI)
        {
            const word& fieldName = fItems_[fieldI];

            const surfaceScalarField& fD =
                obr_.lookupObject<surfaceScalarField>(fieldName);

            const fvMesh& mesh = fD.mesh();

            // 1. integrate over all face zones

            scalarField integralVals(faceZonesSet_.size());

            forAll(faceZonesSet_, setI)
            {
                const word name = faceZonesSet_[setI];

                label zoneID = mesh.faceZones().findZoneID(name);

                const faceZone& fz = mesh.faceZones()[zoneID];

                integralVals[setI] = returnReduce
                (
                    calcFaceZonesIntegral(fD, fz),
                    sumOp<scalar>()
                );
            }


            unsigned int w = IOstream::defaultPrecision() + 7;

            // 2. Write only on master

            if
            (
                Pstream::master()
             && faceZonesIntegrationFilePtr_.found(fieldName)
            )
            {
                OFstream& os = *faceZonesIntegrationFilePtr_(fieldName);

                os  << obr_.time().value();

                forAll(integralVals, setI)
                {
                    os  << ' ' << setw(w) << integralVals[setI];

                    if (log_)
                    {
                        Info<< "faceZonesIntegration output:" << nl
                            << "    Integration" << integralVals[setI] << endl;
                    }
                }

                os  << endl;
            }
        }
    }
}


Foam::scalar Foam::faceZonesIntegration::calcFaceZonesIntegral
(
    const surfaceScalarField& fD,
    const faceZone& fz
) const
{
    scalar dm = 0.0;
    const fvMesh& mesh = fD.mesh();

    forAll (fz, i)
    {
        label faceI = fz[i];

        if (mesh.isInternalFace(faceI))
        {
            if (fz.flipMap()[i])
            {
                dm -= fD[faceI];
            }
            else
            {
                dm += fD[faceI];
            }
        }
        else
        {
            label patchI = mesh.boundaryMesh().whichPatch(faceI);
            const polyPatch& pp = mesh.boundaryMesh()[patchI];
            if (isA<processorPolyPatch>(pp))
            {
                if (refCast<const processorPolyPatch>(pp).owner())
                {
                    if (fz.flipMap()[i])
                    {
                        dm -= fD.boundaryField()[patchI][pp.whichFace(faceI)];
                    }
                    else
                    {
                        dm += fD.boundaryField()[patchI][pp.whichFace(faceI)];
                    }
                }
            }
            else if (isA<cyclicPolyPatch>(pp))
            {
                label patchFaceI = faceI - pp.start();
                if (patchFaceI < pp.size()/2)
                {
                    if (fz.flipMap()[i])
                    {
                        dm -= fD.boundaryField()[patchI][patchFaceI];
                    }
                    else
                    {
                        dm += fD.boundaryField()[patchI][patchFaceI];
                    }
                }
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                label patchFaceI = faceI - pp.start();
                if (fz.flipMap()[i])
                {
                    dm -= fD.boundaryField()[patchI][patchFaceI];
                }
                else
                {
                    dm += fD.boundaryField()[patchI][patchFaceI];
                }
            }
        }
    }

    return dm;
}


// ************************************************************************* //
