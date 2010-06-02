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

#include "vtkFoam.H"

#include "argList.H"
#include "Time.H"
#include "polyBoundaryMeshEntries.H"
#include "IOobjectList.H"
#include "wordList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointMesh.H"
#include "volPointInterpolation.H"

#include "vtkFoamReader.h"
#include "vtkDataArraySelection.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkCharArray.h"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::vtkFoam, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

#include "vtkFoamConvertFields.H"

void Foam::vtkFoam::SetName
(
    vtkUnstructuredGrid* vtkMesh,
    const char* name
)
{
    vtkCharArray* nmArray =  vtkCharArray::New();
    nmArray->SetName("Name");
    size_t len = strlen(name);
    nmArray->SetNumberOfTuples(static_cast<vtkIdType>(len)+1);
    char* copy = nmArray->GetPointer(0);
    memcpy(copy, name, len);
    copy[len] = '\0';
    vtkMesh->GetFieldData()->AddArray(nmArray);
    nmArray->Delete();
}


Foam::string Foam::vtkFoam::padTimeString(const string& ts)
{
    return ts + string("            ", max(label(12 - ts.size()), 0));
}


// Pad the patch name string in order to account for dynamic changes
// in patch names during topological changes
Foam::string Foam::vtkFoam::padPatchString(const string& ps)
{
    label n = max(label(50 - ps.size()), 0);
    return ps + string("                                                  ", n);
}


void Foam::vtkFoam::setSelectedTime
(
    Time& runTime,
    vtkFoamReader* reader
)
{
    // Get times list
    instantList Times = runTime.times();
    int timeIndex = min(max(reader->GetTimeStep() + 1, 0), Times.size()-1);

    // If this is the first call timeIndex will be 0 ("constant")
    // so reset to the first time step if one exists and deselect every
    // element of the selection array
    if (timeIndex == 0)
    {
        timeIndex = min(1, Times.size()-1);
        reader->GetTimeSelection()->DisableAllArrays();
    }

    label selectedTimeIndex = -1;
    label nSelectedTimes = reader->GetTimeSelection()->GetNumberOfArrays();

    for (label i=nSelectedTimes-1; i>=0; i--)
    {
        if(reader->GetTimeSelection()->GetArraySetting(i))
        {
            word timeName = string::validate<word>
            (
                reader->GetTimeSelection()->GetArrayName(i)
            );

            forAll(Times, j)
            {
                if (Times[j].name() == timeName)
                {
                    selectedTimeIndex = j;
                    break;
                }
            }
            break;
        }
    }

    if (selectedTimeIndex != -1)
    {
        timeIndex = min(selectedTimeIndex, Times.size()-1);
    }

    if (debug)
    {
        Info<< "Selecting time " << Times[timeIndex].name() << endl;
    }

    runTime.setTime(Times[timeIndex], timeIndex);

    Times = runTime.times();

    reader->SetTimeStepRange(0, max(Times.size()-2, 0));

    // reset the time steps ...
    reader->GetTimeSelection()->RemoveAllArrays();

    int* TimeStepLimits = reader->GetTimeStepLimits();
    label maxStartTimes = min(Times.size(), TimeStepLimits[0]);
    label maxNTimes = min(Times.size() - maxStartTimes, TimeStepLimits[1]);

    for (label i=0; i<maxStartTimes; i++)
    {
        reader->GetTimeSelection()
            ->AddArray(padTimeString(Times[i].name()).c_str());
    }

    if (Times.size() > TimeStepLimits[0] + TimeStepLimits[1])
    {
        reader->GetTimeSelection()->AddArray(padTimeString("...").c_str());
    }

    for (label i=Times.size() - maxNTimes; i<Times.size(); i++)
    {
        reader->GetTimeSelection()
            ->AddArray(padTimeString(Times[i].name()).c_str());
    }

    // Disable all the time selections (which are all selected by default) ...
    reader->GetTimeSelection()->DisableAllArrays();

    // But maintain the selections made previously
    if (selectedTimeIndex != -1 && selectedTimeIndex < Times.size())
    {
        reader->GetTimeSelection()->EnableArray
            (padTimeString(Times[selectedTimeIndex].name()).c_str());
    }
}


void Foam::vtkFoam::updateSelectedRegions()
{
    if (debug)
    {
        Info<< "Foam::vtkFoam::updateSelectedRegions()" << endl;
    }

    label nRegions = reader_->GetRegionSelection()->GetNumberOfArrays();

    selectedRegions_.setSize(nRegions);

    // Read the selected patches and add to the region list
    for (int i=0; i<nRegions; i++)
    {
        selectedRegions_[i] =
            reader_->GetRegionSelection()->GetArraySetting(i);
    }
}


void Foam::vtkFoam::convertMesh()
{
    if (debug)
    {
        Info<< "Foam::vtkFoam::convertMesh()" << endl;
    }

    const fvMesh& mesh = *meshPtr_;

    // Read the internal mesh as region 0 if selected
    if (reader_->GetRegionSelection()->GetArraySetting(0))
    {
        selectedRegions_[0] = true;
        addInternalMesh
        (
            mesh,
            vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(0))
        );
    }
    else
    {
        selectedRegions_[0] = false;

        vtkUnstructuredGrid *vtkMesh =
            vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(0));

        vtkMesh->Initialize();
        SetName(vtkMesh, "(Internal Mesh)");
    }


    // Read the selected patches and add to the region list

    polyBoundaryMeshEntries patchEntries
    (
        IOobject
        (
            "boundary",
            dbPtr_().findInstance(polyMesh::meshSubDir, "boundary"),
            polyMesh::meshSubDir,
            dbPtr_(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    label regioni = 0;
    label regioniLast = 0;

    // Read in the number Outputs (patch regions) currently being used
    label currNOutputs = reader_->GetNumberOfOutputs();

    // Cycle through all the patches in the boundary file for the relevant
    // time step
    forAll(patchEntries, entryi)
    {
        // Number of faces in the current patch (Used to detect dummy patches
        // of size zero)
        label nFaces(readLabel(patchEntries[entryi].dict().lookup("nFaces")));

        // Check to see if the patch is currently a part of the displayed list
        if
        (
            reader_->GetRegionSelection()->ArrayExists
            (
                padPatchString(patchEntries[entryi].keyword()).c_str()
            )
        )
        {
            if  (!nFaces)
            {
                // Remove patch if it is only a dummy patch in the current
                // time step with zero faces
                reader_->GetRegionSelection()->RemoveArrayByName
                (
                    padPatchString(patchEntries[entryi].keyword()).c_str()
                );
            }
            else
            {
                // A patch already existent in the list and which
                // continues to exist found
                regioni++;
            }
        }
        else
        {
            // A new patch so far not yet included into the list has been found
            if  (nFaces)
            {
                regioni++;

                // Add a new entry to the list of regions
                reader_->GetRegionSelection()->AddArray
                (
                    padPatchString(patchEntries[entryi].keyword()).c_str()
                );

                // AddArray automatically enables a new array... disable
                // it manually
                reader_->GetRegionSelection()->DisableArray
                (
                    padPatchString(patchEntries[entryi].keyword()).c_str()
                );
            }
        }

        // Avoid Initialization of the same Output twice
        if (regioni != regioniLast)
        {
            // Only setup an Output if it has not been setup before
            if(regioni >= currNOutputs)
            {
                vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
                reader_->SetNthOutput(regioni,ugrid);
                ugrid->Delete();
            }
            // Initialize -> Delete memory used, and reset to zero state
            reader_->GetOutput(regioni)->Initialize();
            regioniLast = regioni;
        }
    }

    // Initialize (reset to zero and free) any outputs which are not used
    // anymore
    if (regioni < currNOutputs)
    {
        for(label i = (regioni+1); i < currNOutputs;i++)
        {
            reader_->GetOutput(i)->Initialize();
        }
    }

    selectedRegions_.setSize(regioni + 1);

    regioni = 0;

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll (patches, patchi)
    {
        if (patches[patchi].size())
        {
            regioni++;

            if (reader_->GetRegionSelection()->GetArraySetting(regioni))
            {
                selectedRegions_[regioni] = true;
                addPatch
                (
                    patches[patchi],
                    vtkUnstructuredGrid::SafeDownCast
                    (
                        reader_->GetOutput(regioni)
                    )
                );
            }
            else
            {
                selectedRegions_[regioni] = false;

                vtkUnstructuredGrid *vtkMesh =
                    vtkUnstructuredGrid::SafeDownCast
                    (
                        reader_->GetOutput(regioni)
                    );

                vtkMesh->Initialize();
                SetName
                (
                    vtkMesh,
                    ('(' + padPatchString(patches[patchi].name()) + ')').c_str()
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkFoam::vtkFoam(const char* const FileName, vtkFoamReader* reader)
:
    reader_(reader),
    argsPtr_(NULL),
    dbPtr_(NULL),
    meshPtr_(NULL)
{
    fileName fullCasePath(fileName(FileName).path());

    if (!isDir(fullCasePath))
    {
        return;
    }

    char* argvStrings[3];
    argvStrings[0] = new char[9];
    strcpy(argvStrings[0], "/vtkFoam");
    argvStrings[1] = new char[6];
    strcpy(argvStrings[1], "-case");
    argvStrings[2] = new char[fullCasePath.size()+1];
    strcpy(argvStrings[2], fullCasePath.c_str());

    int argc = 3;
    char** argv = &argvStrings[0];
    argsPtr_.reset(new argList(argc, argv));

    for(int i = 0; i < argc; i++)
    {
        delete[] argvStrings[i];
    }

    dbPtr_.reset
    (
        new Time
        (
            Time::controlDictName,
            argsPtr_().rootPath(),
            argsPtr_().caseName()
        )
    );
    dbPtr_().functionObjects().off();
    setSelectedTime(dbPtr_(), reader_);

    if (debug)
    {
        Info<< "vtkFoam::ExecuteInformation: Initialising outputs" << endl;
    }

    reader_->GetRegionSelection()->AddArray("Internal Mesh");

    vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
    reader_->SetNthOutput(0, ugrid);
    ugrid->Delete();
    reader_->GetOutput(0)->Initialize();

    polyBoundaryMeshEntries patchEntries
    (
        IOobject
        (
            "boundary",
            dbPtr_().findInstance(polyMesh::meshSubDir, "boundary"),
            polyMesh::meshSubDir,
            dbPtr_(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    label regioni = 0;
    forAll(patchEntries, entryi)
    {
        label nFaces(readLabel(patchEntries[entryi].dict().lookup("nFaces")));

        if (nFaces)
        {
            regioni++;

            reader_->GetRegionSelection()->AddArray
            (
                padPatchString(patchEntries[entryi].keyword()).c_str()
            );

            vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
            reader_->SetNthOutput(regioni, ugrid);
            ugrid->Delete();
            reader_->GetOutput(regioni)->Initialize();
        }
    }

    selectedRegions_.setSize(regioni + 1);
    selectedRegions_ = true;

    UpdateInformation();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkFoam::~vtkFoam()
{
    // Do NOT delete meshPtr_ since still referenced somehow.
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "vtkFoamAddFields.H"

void Foam::vtkFoam::UpdateInformation()
{
    if (debug)
    {
        Info<< "TimeStep = " << reader_->GetTimeStep() << endl;
    }

    setSelectedTime(dbPtr_(), reader_);

    // Search for list of objects for this time
    IOobjectList objects(dbPtr_(), dbPtr_().timeName());

    addFields<volScalarField>(reader_->GetVolFieldSelection(), objects);
    addFields<volVectorField>(reader_->GetVolFieldSelection(), objects);
    addFields<volSphericalTensorField>(reader_->GetVolFieldSelection(), objects);
    addFields<volSymmTensorField>(reader_->GetVolFieldSelection(), objects);
    addFields<volTensorField>(reader_->GetVolFieldSelection(), objects);

    addFields<pointScalarField>(reader_->GetPointFieldSelection(), objects);
    addFields<pointVectorField>(reader_->GetPointFieldSelection(), objects);
    addFields<pointSphericalTensorField>(reader_->GetPointFieldSelection(), objects);
    addFields<pointSymmTensorField>(reader_->GetPointFieldSelection(), objects);
    addFields<pointTensorField>(reader_->GetPointFieldSelection(), objects);
}


void Foam::vtkFoam::Update()
{
    if
    (
        !reader_->GetCacheMesh()
     || reader_->GetTimeSelection()->GetArraySetting(0)
    )
    {
        meshPtr_= NULL;
    }

    // Clear the current set of selected fields

    for (label i=0; i<reader_->GetNumberOfOutputs(); i++)
    {
        vtkUnstructuredGrid *vtkMesh =
            vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(i));

        vtkCellData* cellData = vtkMesh->GetCellData();
        int numberOfCellArrays = cellData->GetNumberOfArrays();

        wordList cellFieldNames(numberOfCellArrays);
        for (int j=0; j<numberOfCellArrays; j++)
        {
            cellFieldNames[j] = cellData->GetArrayName(j);
        }

        for (int j=0; j<numberOfCellArrays; j++)
        {
            cellData->RemoveArray(cellFieldNames[j].c_str());
        }

        vtkPointData* pointData = vtkMesh->GetPointData();
        int numberOfPointArrays = pointData->GetNumberOfArrays();

        wordList pointFieldNames(numberOfPointArrays);
        for (int j=0; j<numberOfPointArrays; j++)
        {
            pointFieldNames[j] = pointData->GetArrayName(j);
        }

        for (int j=0; j<numberOfPointArrays; j++)
        {
            pointData->RemoveArray(pointFieldNames[j].c_str());
        }
    }

    // Check to see if the mesh has been created

    if (!meshPtr_)
    {
        if (debug)
        {
            Info<< "Reading Mesh" << endl;
        }
        meshPtr_ =
            new fvMesh
            (
                IOobject
                (
                    fvMesh::defaultRegion,
                    dbPtr_().timeName(),
                    dbPtr_()
                )
            );
        convertMesh();
    }
    else
    {
        boolList oldSelectedRegions = selectedRegions_;
        updateSelectedRegions();
        if
        (
            meshPtr_->readUpdate() != fvMesh::UNCHANGED
         || oldSelectedRegions != selectedRegions_
        )
        {
            convertMesh();
        }
    }

    if (debug)
    {
        Info<< "converting fields" << endl;
    }

    const fvMesh& mesh = *meshPtr_;

    Foam::volPointInterpolation pInterp(mesh);

    // Search for list of objects for this time
    Foam::IOobjectList objects(mesh, dbPtr_().timeName());

    convertVolFields<Foam::scalar>
    (
        mesh, pInterp, objects, reader_->GetVolFieldSelection()
    );
    convertVolFields<Foam::vector>
    (
        mesh, pInterp, objects, reader_->GetVolFieldSelection()
    );
    convertVolFields<Foam::sphericalTensor>
    (
        mesh, pInterp, objects, reader_->GetVolFieldSelection()
    );
    convertVolFields<Foam::symmTensor>
    (
        mesh, pInterp, objects, reader_->GetVolFieldSelection()
    );
    convertVolFields<Foam::tensor>
    (
        mesh, pInterp, objects, reader_->GetVolFieldSelection()
    );

    convertPointFields<Foam::scalar>
    (
        mesh, objects, reader_->GetPointFieldSelection()
    );
    convertPointFields<Foam::vector>
    (
        mesh, objects, reader_->GetPointFieldSelection()
    );
    convertPointFields<Foam::sphericalTensor>
    (
        mesh, objects, reader_->GetPointFieldSelection()
    );
    convertPointFields<Foam::symmTensor>
    (
        mesh, objects, reader_->GetPointFieldSelection()
    );
    convertPointFields<Foam::tensor>
    (
        mesh, objects, reader_->GetPointFieldSelection()
    );

    if (debug)
    {
        Info<< "done" << endl;
    }
}


// ************************************************************************* //
