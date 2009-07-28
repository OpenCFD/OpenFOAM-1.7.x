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

InClass
    vtkPV3Foam

\*---------------------------------------------------------------------------*/

#ifndef vtkPV3FoamFaceField_H
#define vtkPV3FoamFaceField_H

// VTK includes
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkPolyData.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkPV3Foam::convertFaceField
(
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const partInfo& selector,
    const label datasetNo,
    const fvMesh& mesh,
    const labelList& faceLabels
)
{
    const label nComp = pTraits<Type>::nComponents;
    const label nInternalFaces = mesh.nInternalFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();

    vtkFloatArray *cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples( faceLabels.size() );
    cellData->SetNumberOfComponents( nComp );
    cellData->Allocate( nComp*faceLabels.size() );
    cellData->SetName( tf.name().c_str() );

    if (debug)
    {
        Info<< "convert convertFaceField: "
            << tf.name()
            << " size = " << tf.size()
            << " nComp=" << nComp
            << " nTuples = " << faceLabels.size() <<  endl;
    }

    float vec[nComp];

    // for interior faces: average owner/neighbour
    // for boundary faces: owner
    forAll(faceLabels, faceI)
    {
        const label faceNo = faceLabels[faceI];
        if (faceNo < nInternalFaces)
        {
            Type t = 0.5*(tf[faceOwner[faceNo]] + tf[faceNeigh[faceNo]]);

            for (direction d=0; d<nComp; d++)
            {
                vec[d] = component(t, d);
            }
        }
        else
        {
            const Type& t = tf[faceOwner[faceNo]];
            for (direction d=0; d<nComp; d++)
            {
                vec[d] = component(t, d);
            }
        }

        cellData->InsertTuple(faceI, vec);
    }


    vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    )   ->GetCellData()
        ->AddArray(cellData);

    cellData->Delete();
}


template<class Type>
void Foam::vtkPV3Foam::convertFaceField
(
    const GeometricField<Type, fvPatchField, volMesh>& tf,
    vtkMultiBlockDataSet* output,
    const partInfo& selector,
    const label datasetNo,
    const fvMesh& mesh,
    const faceSet& fSet
)
{
    const label nComp = pTraits<Type>::nComponents;
    const label nInternalFaces = mesh.nInternalFaces();
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();

    vtkFloatArray *cellData = vtkFloatArray::New();
    cellData->SetNumberOfTuples( fSet.size() );
    cellData->SetNumberOfComponents( nComp );
    cellData->Allocate( nComp*fSet.size() );
    cellData->SetName( tf.name().c_str() );

    if (debug)
    {
        Info<< "convert convertFaceField: "
            << tf.name()
            << " size = " << tf.size()
            << " nComp=" << nComp
            << " nTuples = " << fSet.size() <<  endl;
    }

    float vec[nComp];

    // for interior faces: average owner/neighbour
    // for boundary faces: owner
    label faceI = 0;
    forAllConstIter(faceSet, fSet, iter)
    {
        const label faceNo = iter.key();

        if (faceNo < nInternalFaces)
        {
            Type t = 0.5*(tf[faceOwner[faceNo]] + tf[faceNeigh[faceNo]]);

            for (direction d=0; d<nComp; d++)
            {
                vec[d] = component(t, d);
            }
        }
        else
        {
            const Type& t = tf[faceOwner[faceNo]];
            for (direction d=0; d<nComp; d++)
            {
                vec[d] = component(t, d);
            }
        }

        cellData->InsertTuple(faceI, vec);
        ++faceI;
    }


    vtkPolyData::SafeDownCast
    (
        GetDataSetFromBlock(output, selector, datasetNo)
    )   ->GetCellData()
        ->AddArray(cellData);

    cellData->Delete();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
