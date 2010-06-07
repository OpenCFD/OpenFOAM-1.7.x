/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include <ctype.h>

#include "vtkFoamReader.h"

#include "vtkCallbackCommand.h"
#include "vtkDataArraySelection.h"
#include "vtkDataArrayCollection.h"
#include "vtkObjectFactory.h"
#include "vtkDataSet.h"
#include "vtkErrorCode.h"
#include "vtkUnstructuredGrid.h"

#include "vtkFoam.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

vtkCxxRevisionMacro(vtkFoamReader, "$Revision: 1.20 $");
vtkStandardNewMacro(vtkFoamReader);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vtkFoamReader::vtkFoamReader()
{
    StoredOutputs = NULL;

    FileName  = NULL;
    foamData_ = NULL;

    CacheMesh = 0;

    UpdateGUI = 1;
    UpdateGUIOld = 1;
    TimeStep = 0;
    TimeStepRange[0] = 0;
    TimeStepRange[1] = 0;

    TimeStepLimits[0] = 2;
    TimeStepLimits[1] = 5;

    TimeSelection = vtkDataArraySelection::New();
    RegionSelection = vtkDataArraySelection::New();
    VolFieldSelection = vtkDataArraySelection::New();
    PointFieldSelection = vtkDataArraySelection::New();

    // Setup the selection callback to modify this object when an array
    // selection is changed.
    SelectionObserver = vtkCallbackCommand::New();
    SelectionObserver->SetCallback(&vtkFoamReader::SelectionModifiedCallback);
    SelectionObserver->SetClientData(this);

    TimeSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
    RegionSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
    VolFieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );
    PointFieldSelection->AddObserver
    (
        vtkCommand::ModifiedEvent,
        this->SelectionObserver
    );

    // This is needed by ParaView 2.?.?
    this->SetNumberOfOutputPorts(0);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

vtkFoamReader::~vtkFoamReader()
{
    if (foamData_)
    {
        delete foamData_;
    }

    if (StoredOutputs)
    {
        StoredOutputs->Delete();
    }

    if (FileName)
    {
        delete [] FileName;
    }

    TimeSelection->RemoveObserver(this->SelectionObserver);
    RegionSelection->RemoveObserver(this->SelectionObserver);
    VolFieldSelection->RemoveObserver(this->SelectionObserver);
    PointFieldSelection->RemoveObserver(this->SelectionObserver);
    SelectionObserver->Delete();

    TimeSelection->Delete();
    RegionSelection->Delete();
    VolFieldSelection->Delete();
    PointFieldSelection->Delete();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void vtkFoamReader::ExecuteInformation()
{
    if (!foamData_)
    {
        vtkDebugMacro( << "Reading Foam case" << FileName);
        foamData_ = new Foam::vtkFoam(FileName, this);
    }
    else
    {
        foamData_->UpdateInformation();
    }

    vtkDebugMacro( << "end of ExecuteInformation\n");
}


void vtkFoamReader::Execute()
{
    if (!StoredOutputs)
    {
        foamData_->Update();

        StoredOutputs = vtkFoamData::New();

        for (int i = 0; i < GetNumberOfOutputs(); i++)
        {
            vtkDataObject* tmp = GetOutput(i);
            vtkDataObject* output = tmp->NewInstance();
            output->ShallowCopy(tmp);
            StoredOutputs->SetNthOutput(i, output);
            output->Delete();
        }
    }
    else
    {
        for (int i = 0; i < GetNumberOfOutputs(); i++)
        {
            vtkDataObject* output = GetOutput(i);
            int tempExtent[6];
            output->GetUpdateExtent(tempExtent);
            output->ShallowCopy(StoredOutputs->GetOutput(i));
            output->SetUpdateExtent(tempExtent);
        }

        if (UpdateGUIOld == GetUpdateGUI())
        {
            foamData_->Update();

            for (int i = 0; i < GetNumberOfOutputs(); i++)
            {
                vtkDataObject* tmp = GetOutput(i);
                vtkDataObject* output = tmp->NewInstance();
                output->ShallowCopy(tmp);
                StoredOutputs->SetNthOutput(i, output);
                output->Delete();
            }
        }
    }

    UpdateGUIOld = GetUpdateGUI();
}


void vtkFoamReader::SetFileName(const char *name)
{
    if (name && !FileName || (FileName && !strcmp(FileName,name)))
    {
        if (!FileName)
        {
            FileName = new char[strlen(name) + 1];
            strcpy(FileName, name);
        }
    }
    else
    {
        vtkErrorMacro("Changing case is not currently supported.\nPlease delete reader and create a new one for the new case.");
        return;
    }

    /*
    if ( FileName && name && (!strcmp(FileName,name)))
    {
        return;
    }

    if (!name && !FileName)
    {
        return;
    }

    if (FileName)
    {
        delete [] FileName;
    }

    FileName = new char[strlen(name) + 1];
    strcpy(FileName, name);

    if (foamData_)
    {
        delete foamData_;
        foamData_ = NULL;

        if (StoredOutputs)
        {
            StoredOutputs->Delete();
            StoredOutputs = NULL;
        }
    }

    Modified();
    */
}


void vtkFoamReader::PrintSelf(ostream& os, vtkIndent indent)
{
    Superclass::PrintSelf(os,indent);

    os  << indent << "File Name: " 
        << (FileName ? FileName : "(none)") << "\n";
}


vtkDataArraySelection* vtkFoamReader::GetTimeSelection()
{
    return TimeSelection;
}

int vtkFoamReader::GetNumberOfTimeArrays()
{
    return TimeSelection->GetNumberOfArrays();
}

const char* vtkFoamReader::GetTimeArrayName(int index)
{
    return TimeSelection->GetArrayName(index);
}

int vtkFoamReader::GetTimeArrayStatus(const char* name)
{
    return TimeSelection->ArrayIsEnabled(name);
}

void vtkFoamReader::SetTimeArrayStatus(const char* name, int status)
{
    if(status)
    {
        TimeSelection->EnableArray(name);
    }
    else
    {
        TimeSelection->DisableArray(name);
    }
}

vtkDataArraySelection* vtkFoamReader::GetRegionSelection()
{
    return RegionSelection;
}

int vtkFoamReader::GetNumberOfRegionArrays()
{
    return RegionSelection->GetNumberOfArrays();
}

const char* vtkFoamReader::GetRegionArrayName(int index)
{
    return RegionSelection->GetArrayName(index);
}

int vtkFoamReader::GetRegionArrayStatus(const char* name)
{
    return RegionSelection->ArrayIsEnabled(name);
}

void vtkFoamReader::SetRegionArrayStatus(const char* name, int status)
{
    if(status)
    {
        RegionSelection->EnableArray(name);
    }
    else
    {
        RegionSelection->DisableArray(name);
    }
}


vtkDataArraySelection* vtkFoamReader::GetVolFieldSelection()
{
    return VolFieldSelection;
}

int vtkFoamReader::GetNumberOfVolFieldArrays()
{
    return VolFieldSelection->GetNumberOfArrays();
}

const char* vtkFoamReader::GetVolFieldArrayName(int index)
{
    return VolFieldSelection->GetArrayName(index);
}

int vtkFoamReader::GetVolFieldArrayStatus(const char* name)
{
    return VolFieldSelection->ArrayIsEnabled(name);
}

void vtkFoamReader::SetVolFieldArrayStatus(const char* name, int status)
{
    if(status)
    {
        VolFieldSelection->EnableArray(name);
    }
    else
    {
        VolFieldSelection->DisableArray(name);
    }
}


vtkDataArraySelection* vtkFoamReader::GetPointFieldSelection()
{
    return PointFieldSelection;
}

int vtkFoamReader::GetNumberOfPointFieldArrays()
{
    return PointFieldSelection->GetNumberOfArrays();
}

const char* vtkFoamReader::GetPointFieldArrayName(int index)
{
    return PointFieldSelection->GetArrayName(index);
}

int vtkFoamReader::GetPointFieldArrayStatus(const char* name)
{
    return PointFieldSelection->ArrayIsEnabled(name);
}

void vtkFoamReader::SetPointFieldArrayStatus(const char* name, int status)
{
    if(status)
    {
        PointFieldSelection->EnableArray(name);
    }
    else
    {
        PointFieldSelection->DisableArray(name);
    }
}


void vtkFoamReader::SelectionModifiedCallback
(
    vtkObject*,
    unsigned long,
    void* clientdata,
    void*
)
{
    static_cast<vtkFoamReader*>(clientdata)->SelectionModified();
}

void vtkFoamReader::SelectionModified()
{
    Modified();
}


// ************************************************************************* //
