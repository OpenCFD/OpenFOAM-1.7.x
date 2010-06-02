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

Class
    vtkFoamReader

Description

SourceFiles
    vtkFoamReader.cxx

\*---------------------------------------------------------------------------*/

#ifndef vtkFoamReader_h
#define vtkFoamReader_h

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "vtkDataSetSource.h"
#include "vtkFoamData.h"

// * * * * * * * * * * * * * Forward Declarations  * * * * * * * * * * * * * //

namespace Foam
{
    class vtkFoam;
}

class vtkPoints;
class vtkDataArraySelection;
class vtkDataArrayCollection;
class vtkCallbackCommand;

/*---------------------------------------------------------------------------*\
                         Class vtkFoamReader Declaration
\*---------------------------------------------------------------------------*/

class VTK_IO_EXPORT vtkFoamReader
: 
    public vtkDataSetSource
{

public:

    //- Standard VTK class creation function
    static vtkFoamReader *New();

    //- Standard VTK class type and revision declaration macro
    vtkTypeRevisionMacro(vtkFoamReader,vtkDataSetSource);

    //- Standard VTK class print function
    void PrintSelf(ostream& os, vtkIndent indent);

    // File name of FOAM datafile to read
    void SetFileName(const char *);
    //vtkSetStringMacro(FileName);
    vtkGetStringMacro(FileName);

    // GUI update control
    vtkSetMacro(UpdateGUI, int);
    vtkGetMacro(UpdateGUI, int);

    // FOAM mesh caching control
    vtkSetMacro(CacheMesh, int);
    vtkGetMacro(CacheMesh, int);

    // Time-step slider control
    vtkSetMacro(TimeStep, int);
    vtkGetMacro(TimeStep, int);
    vtkSetVector2Macro(TimeStepRange, int);
    vtkGetVector2Macro(TimeStepRange, int);

    // Control of the upper and lower limits on the number of times
    // displayed in the selection list
    vtkSetVector2Macro(TimeStepLimits, int);
    vtkGetVector2Macro(TimeStepLimits, int);

    // Time selection list control
    vtkDataArraySelection* GetTimeSelection();
    int GetNumberOfTimeArrays();
    const char* GetTimeArrayName(int index);
    int GetTimeArrayStatus(const char* name);
    void SetTimeArrayStatus(const char* name, int status);

    // Region selection list control
    vtkDataArraySelection* GetRegionSelection();
    int GetNumberOfRegionArrays();
    const char* GetRegionArrayName(int index);
    int GetRegionArrayStatus(const char* name);
    void SetRegionArrayStatus(const char* name, int status);

    // volField selection list control
    vtkDataArraySelection* GetVolFieldSelection();
    int GetNumberOfVolFieldArrays();
    const char* GetVolFieldArrayName(int index);
    int GetVolFieldArrayStatus(const char* name);
    void SetVolFieldArrayStatus(const char* name, int status);

    // pointField selection list control
    vtkDataArraySelection* GetPointFieldSelection();
    int GetNumberOfPointFieldArrays();
    const char* GetPointFieldArrayName(int index);
    int GetPointFieldArrayStatus(const char* name);
    void SetPointFieldArrayStatus(const char* name, int status);  

    // SetNthOutput provided so that vtkFoam can access it
    void SetNthOutput(int num, vtkDataObject *output)
    {
        vtkDataSetSource::SetNthOutput(num, output);
    }

    // Standard VTK ExecuteInformation function overriding the base-class.
    // Called by ParaView before GUI is displayed.
    virtual void ExecuteInformation();

    // Callback registered with the SelectionObserver
    // for all the selection lists
    static void SelectionModifiedCallback
    (
        vtkObject* caller,
        unsigned long eid,
        void* clientdata,
        void* calldata
    );

    void SelectionModified();


protected:

    vtkFoamReader();
    ~vtkFoamReader();

    // Standard VTK execute function overriding the base-class.
    // Called by ParaView when Accept is pressed.
    void Execute();

    // Cache for the outputs.  These are stored before the end of Execute()
    // and re-instated at the beginning because the Outputs would disappear
    // otherwise.
    vtkFoamData* StoredOutputs;

    // FOAM file name (*.foam)
    char *FileName;

    //BTX
    Foam::vtkFoam* foamData_;
    //ETX

    int CacheMesh;

    int UpdateGUI;
    int UpdateGUIOld;
    int TimeStep;
    int TimeStepRange[2];

    int TimeStepLimits[2];

    vtkDataArraySelection* TimeSelection;
    vtkDataArraySelection* RegionSelection;
    vtkDataArraySelection* VolFieldSelection;
    vtkDataArraySelection* PointFieldSelection;

    // The observer to modify this object when the array selections are modified
    vtkCallbackCommand* SelectionObserver;


private:

    vtkFoamReader(const vtkFoamReader&);  // Not implemented.
    void operator=(const vtkFoamReader&);  // Not implemented.
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
