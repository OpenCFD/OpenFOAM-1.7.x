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

\*---------------------------------------------------------------------------*/

#include "vtkPV3Foam.H"

// Foam includes
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "polyBoundaryMeshEntries.H"
#include "entry.H"
#include "Cloud.H"
#include "vtkPV3FoamReader.h"

// local headers
#include "vtkPV3FoamAddToSelection.H"
#include "vtkPV3FoamUpdateInfoFields.H"

// VTK includes
#include "vtkDataArraySelection.h"


// * * * * * * * * * * * * * * * Private Classes * * * * * * * * * * * * * * //

namespace Foam
{

//- A class for reading zone information without requiring a mesh
class zonesEntries
:
    public regIOobject,
    public PtrList<entry>
{

public:

    // Constructors

        explicit zonesEntries(const IOobject& io)
        :
            regIOobject(io),
            PtrList<entry>(readStream("regIOobject"))
        {
            close();
        }

   // Member functions

        bool writeData(Ostream&) const
        {
            notImplemented("zonesEntries::writeData(Ostream&) const");
            return false;
        }
};

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList Foam::vtkPV3Foam::readZoneNames(const word& zoneType)
{
    wordList zoneNames;

    // mesh not loaded - read from file
    IOobject ioObj
    (
        zoneType,
        dbPtr_().findInstance
        (
            meshDir_,
            zoneType,
            IOobject::READ_IF_PRESENT
        ),
        meshDir_,
        dbPtr_(),
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (ioObj.headerOk())
    {
        zonesEntries zones(ioObj);

        zoneNames.setSize(zones.size());
        forAll(zones, zoneI)
        {
            zoneNames[zoneI] = zones[zoneI].keyword();
        }
    }

    return zoneNames;
}


void Foam::vtkPV3Foam::updateInfoInternalMesh()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoInternalMesh" << endl;
    }

    vtkDataArraySelection* partSelection = reader_->GetPartSelection();

    // Determine mesh parts (internalMesh, patches...)
    //- Add internal mesh as first entry
    partInfoVolume_ = partSelection->GetNumberOfArrays();
    partSelection->AddArray("internalMesh");
    partInfoVolume_ += 1;

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(partSelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoInternalMesh" << endl;
    }

}


void Foam::vtkPV3Foam::updateInfoLagrangian()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoLagrangian" << nl
            << "    " << dbPtr_->timePath()/cloud::prefix << endl;
    }


    // use the db directly since this might be called without a mesh,
    // but the region must get added back in
    fileName lagrangianPrefix(cloud::prefix);
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        lagrangianPrefix = meshRegion_/cloud::prefix;
    }

    // Search for list of lagrangian objects for this time
    fileNameList cloudDirs
    (
        readDir(dbPtr_->timePath()/lagrangianPrefix, fileName::DIRECTORY)
    );

    vtkDataArraySelection* partSelection = reader_->GetPartSelection();
    partInfoLagrangian_ = partSelection->GetNumberOfArrays();

    int nClouds = 0;
    forAll(cloudDirs, cloudI)
    {
        // Add cloud to GUI list
        partSelection->AddArray
        (
            (cloudDirs[cloudI] + " - lagrangian").c_str()
        );

        ++nClouds;
    }

    partInfoLagrangian_ += nClouds;

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(partSelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoLagrangian" << endl;
    }
}


void Foam::vtkPV3Foam::updateInfoPatches()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoPatches"
            << " [meshPtr=" << (meshPtr_ ? "set" : "NULL") << "]" << endl;
    }

    vtkDataArraySelection* partSelection = reader_->GetPartSelection();
    partInfoPatches_ = partSelection->GetNumberOfArrays();

    int nPatches = 0;
    if (meshPtr_)
    {
        const polyBoundaryMesh& patches = meshPtr_->boundaryMesh();
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.size())
            {
                // Add patch to GUI list
                partSelection->AddArray
                (
                    (pp.name() + " - patch").c_str()
                );

                ++nPatches;
            }
        }
    }
    else
    {
        // mesh not loaded - read from file
        // but this could fail if we've supplied a bad region name
        IOobject ioObj
        (
            "boundary",
            dbPtr_().findInstance
            (
                meshDir_,
                "boundary",
                IOobject::READ_IF_PRESENT
            ),
            meshDir_,
            dbPtr_(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        // this should only ever fail if the mesh region doesn't exist
        if (ioObj.headerOk())
        {
            polyBoundaryMeshEntries patchEntries(ioObj);

            // Add (non-zero) patches to the list of mesh parts
            forAll(patchEntries, entryI)
            {
                label nFaces
                (
                    readLabel(patchEntries[entryI].dict().lookup("nFaces"))
                );

                // Valid patch if nFace > 0 - add patch to GUI list
                if (nFaces)
                {
                    partSelection->AddArray
                    (
                        (patchEntries[entryI].keyword() + " - patch").c_str()
                    );

                    ++nPatches;
                }
            }
        }
    }
    partInfoPatches_ += nPatches;

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(partSelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoPatches" << endl;
    }
}


void Foam::vtkPV3Foam::updateInfoZones()
{
    if (!reader_->GetIncludeZones())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoZones"
            << " [meshPtr=" << (meshPtr_ ? "set" : "NULL") << "]" << endl;
    }

    vtkDataArraySelection* partSelection = reader_->GetPartSelection();
    wordList namesLst;

    //
    // cellZones information
    // ~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = meshPtr_->cellZones().names();
    }
    else
    {
        namesLst = readZoneNames("cellZones");
    }

    partInfoCellZones_ = partSelection->GetNumberOfArrays();
    forAll(namesLst, elemI)
    {
        partSelection->AddArray((namesLst[elemI] + " - cellZone").c_str());
    }
    partInfoCellZones_ += namesLst.size();


    //
    // faceZones information
    // ~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = meshPtr_->faceZones().names();
    }
    else
    {
        namesLst = readZoneNames("faceZones");
    }

    partInfoFaceZones_ = partSelection->GetNumberOfArrays();
    forAll(namesLst, elemI)
    {
        partSelection->AddArray
        (
            (namesLst[elemI] + " - faceZone").c_str()
        );
    }
    partInfoFaceZones_ += namesLst.size();


    //
    // pointZones information
    // ~~~~~~~~~~~~~~~~~~~~~~
    if (meshPtr_)
    {
        namesLst = meshPtr_->pointZones().names();
    }
    else
    {
        namesLst = readZoneNames("pointZones");
    }

    partInfoPointZones_ = partSelection->GetNumberOfArrays();
    forAll(namesLst, elemI)
    {
        partSelection->AddArray
        (
            (namesLst[elemI] + " - pointZone").c_str()
        );
    }
    partInfoPointZones_ += namesLst.size();


    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(partSelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoZones" << endl;
    }
}


void Foam::vtkPV3Foam::updateInfoSets()
{
    if (!reader_->GetIncludeSets())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoSets" << endl;
    }

    vtkDataArraySelection* partSelection = reader_->GetPartSelection();

    // Add names of sets
    IOobjectList objects
    (
        dbPtr_(),
        dbPtr_().findInstance(meshDir_, "faces", IOobject::READ_IF_PRESENT),
        meshDir_/"sets"
    );


    partInfoCellSets_ = partSelection->GetNumberOfArrays();
    partInfoCellSets_ += addToSelection<cellSet>
    (
        partSelection,
        objects,
        " - cellSet"
    );

    partInfoFaceSets_ = partSelection->GetNumberOfArrays();
    partInfoFaceSets_ += addToSelection<faceSet>
    (
        partSelection,
        objects,
        " - faceSet"
    );

    partInfoPointSets_ = partSelection->GetNumberOfArrays();
    partInfoPointSets_ += addToSelection<pointSet>
    (
        partSelection,
        objects,
        " - pointSet"
    );

    if (debug)
    {
        // just for debug info
        getSelectedArrayEntries(partSelection);

        Info<< "<end> Foam::vtkPV3Foam::updateInfoSets" << endl;
    }
}


void Foam::vtkPV3Foam::updateInfoLagrangianFields()
{
    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::updateInfoLagrangianFields"
            << endl;
    }

    vtkDataArraySelection *fieldSelection =
        reader_->GetLagrangianFieldSelection();

    // preserve the enabled selections
    stringList enabledEntries = getSelectedArrayEntries(fieldSelection);
    fieldSelection->RemoveAllArrays();

    //
    // TODO - currently only get fields from ONE cloud
    // have to decide if the second set of fields get mixed in
    // or dealt with separately

    const partInfo& selector = partInfoLagrangian_;
    int partId = selector.start();

    if (!selector.size() || partId < 0)
    {
        return;
    }

    word cloudName = getPartName(partId);

    // use the db directly since this might be called without a mesh,
    // but the region must get added back in
    fileName lagrangianPrefix(cloud::prefix);
    if (meshRegion_ != polyMesh::defaultRegion)
    {
        lagrangianPrefix = meshRegion_/cloud::prefix;
    }

    IOobjectList objects
    (
        dbPtr_(),
        dbPtr_().timeName(),
        lagrangianPrefix/cloudName
    );

    addToSelection<IOField<label> >
    (
        fieldSelection,
        objects
    );
    addToSelection<IOField<scalar> >
    (
        fieldSelection,
        objects
    );
    addToSelection<IOField<vector> >
    (
        fieldSelection,
        objects
    );
    addToSelection<IOField<sphericalTensor> >
    (
        fieldSelection,

        objects
    );
    addToSelection<IOField<symmTensor> >
    (
        fieldSelection,
        objects
    );
    addToSelection<IOField<tensor> >
    (
        fieldSelection,
        objects
    );

    // restore the enabled selections
    setSelectedArrayEntries(fieldSelection, enabledEntries);

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::updateInfoLagrangianFields - "
            << "lagrangian objects.size() = " << objects.size() << endl;
    }
}


// ************************************************************************* //
