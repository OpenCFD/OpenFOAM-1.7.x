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

Description

\*---------------------------------------------------------------------------*/

#include "vtkPV3Foam.H"

// Foam includes
#include "IOobjectList.H"
#include "vtkPV3FoamReader.h"

// VTK includes
#include "vtkDataArraySelection.h"
#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

#include "vtkPV3FoamVolFields.H"
#include "vtkPV3FoamPointFields.H"
#include "vtkPV3FoamLagrangianFields.H"


void Foam::vtkPV3Foam::pruneObjectList
(
    IOobjectList& objects,
    const wordHashSet& selected
)
{
    // hash all the selected field names
    if (!selected.size())
    {
        objects.clear();
    }

    // only keep selected fields
    forAllIter(IOobjectList, objects, iter)
    {
        if (!selected.found(iter()->name()))
        {
            objects.erase(iter);
        }
    }
}


void Foam::vtkPV3Foam::convertVolFields
(
    vtkMultiBlockDataSet* output
)
{
    const fvMesh& mesh = *meshPtr_;

    wordHashSet selectedFields = getSelected
    (
        reader_->GetVolFieldSelection()
    );

    if (!selectedFields.size())
    {
        return;
    }

    // Get objects (fields) for this time - only keep selected fields
    // the region name is already in the mesh db
    IOobjectList objects(mesh, dbPtr_().timeName());
    pruneObjectList(objects, selectedFields);

    if (!objects.size())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertVolFields" << nl
            << "converting Foam volume fields" << endl;
        forAllConstIter(IOobjectList, objects, iter)
        {
            Info<< "  " << iter()->name()
                << " == " << iter()->objectPath() << nl;
        }
        printMemory();
    }


    PtrList<PrimitivePatchInterpolation<primitivePatch> >
        ppInterpList(mesh.boundaryMesh().size());

    forAll(ppInterpList, i)
    {
        ppInterpList.set
        (
            i,
            new PrimitivePatchInterpolation<primitivePatch>
            (
                mesh.boundaryMesh()[i]
            )
        );
    }


    convertVolFields<scalar>
    (
        mesh, ppInterpList, objects, output
    );
    convertVolFields<vector>
    (
        mesh, ppInterpList, objects, output
    );
    convertVolFields<sphericalTensor>
    (
        mesh, ppInterpList, objects, output
    );
    convertVolFields<symmTensor>
    (
        mesh, ppInterpList, objects, output
    );
    convertVolFields<tensor>
    (
        mesh, ppInterpList, objects, output
    );

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertVolFields" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertPointFields
(
    vtkMultiBlockDataSet* output
)
{
    const fvMesh& mesh = *meshPtr_;

    wordHashSet selectedFields = getSelected
    (
        reader_->GetPointFieldSelection()
    );

    if (!selectedFields.size())
    {
        return;
    }

    // Get objects (fields) for this time - only keep selected fields
    // the region name is already in the mesh db
    IOobjectList objects(mesh, dbPtr_().timeName());
    pruneObjectList(objects, selectedFields);

    if (!objects.size())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertPointFields" << nl
            << "converting Foam volume fields" << endl;
        forAllConstIter(IOobjectList, objects, iter)
        {
            Info<< "  " << iter()->name()
                << " == " << iter()->objectPath() << nl;
        }
        printMemory();
    }

    // Construct interpolation on the raw mesh
    pointMesh pMesh(mesh);


    convertPointFields<scalar>
    (
        mesh, pMesh, objects, output
    );
    convertPointFields<vector>
    (
        mesh, pMesh, objects, output
    );
    convertPointFields<sphericalTensor>
    (
        mesh, pMesh, objects, output
    );
    convertPointFields<symmTensor>
    (
        mesh, pMesh, objects, output
    );
    convertPointFields<tensor>
    (
        mesh, pMesh, objects, output
    );

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertPointFields" << endl;
        printMemory();
    }
}


void Foam::vtkPV3Foam::convertLagrangianFields
(
    vtkMultiBlockDataSet* output
)
{
    partInfo& selector = partInfoLagrangian_;
    const fvMesh& mesh = *meshPtr_;

    wordHashSet selectedFields = getSelected
    (
        reader_->GetLagrangianFieldSelection()
    );

    if (!selectedFields.size())
    {
        return;
    }

    if (debug)
    {
        Info<< "<beg> Foam::vtkPV3Foam::convertLagrangianFields" << endl;
        printMemory();
    }

    for (int partId = selector.start(); partId < selector.end(); ++partId)
    {
        const word  cloudName = getPartName(partId);
        const label datasetNo = partDataset_[partId];

        if (!partStatus_[partId] || datasetNo < 0)
        {
            continue;
        }


        // Get the Lagrangian fields for this time and this cloud
        // but only keep selected fields
        // the region name is already in the mesh db
        IOobjectList objects
        (
            mesh,
            dbPtr_().timeName(),
            cloud::prefix/cloudName
        );
        pruneObjectList(objects, selectedFields);

        if (!objects.size())
        {
            continue;
        }

        if (debug)
        {
            Info<< "converting Foam lagrangian fields" << nl;
            forAllConstIter(IOobjectList, objects, iter)
            {
                Info<< "  " << iter()->name()
                    << " == " << iter()->objectPath() << nl;
            }
        }

        convertLagrangianFields<label>
        (
            objects, output, datasetNo
        );
        convertLagrangianFields<scalar>
        (
            objects, output, datasetNo
        );
        convertLagrangianFields<vector>
        (
            objects, output, datasetNo
        );
        convertLagrangianFields<sphericalTensor>
        (
            objects, output, datasetNo
        );
        convertLagrangianFields<symmTensor>
        (
            objects, output, datasetNo
        );
        convertLagrangianFields<tensor>
        (
            objects, output, datasetNo
        );
    }

    if (debug)
    {
        Info<< "<end> Foam::vtkPV3Foam::convertLagrangianFields" << endl;
        printMemory();
    }
}


// ************************************************************************* //
