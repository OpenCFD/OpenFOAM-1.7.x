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
    Foam::vtkFoam

\*---------------------------------------------------------------------------*/

#ifndef vtkFoamConvertVolField_H
#define vtkFoamConvertVolField_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkFoam::convertVolField
(
    const GeometricField<Type, fvPatchField, volMesh>& tf
)
{
    vtkUnstructuredGrid *vtkMesh =
        vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(0));

    vtkFloatArray *cellTypes = vtkFloatArray::New();
    cellTypes->SetNumberOfTuples(superCells_.size());
    cellTypes->SetNumberOfComponents(Type::nComponents);
    cellTypes->Allocate(Type::nComponents*tf.size());
    cellTypes->SetName(tf.name().c_str());

    float vec[Type::nComponents];

    forAll(superCells_, sci)
    {
        const Type& t = tf[superCells_[sci]];
        for (direction d=0; d<Type::nComponents; d++)
        {
            vec[d] = t[d];
        }

        cellTypes->InsertTuple(sci, vec);
    }

    vtkMesh->GetCellData()->AddArray(cellTypes);
    cellTypes->Delete();
}


template<>
void Foam::vtkFoam::convertVolField
(
    const GeometricField<scalar, fvPatchField, volMesh>& sf
)
{
    vtkUnstructuredGrid *vtkMesh =
        vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(0));

    vtkFloatArray *cellScalars = vtkFloatArray::New();
    cellScalars->SetNumberOfTuples(superCells_.size());
    cellScalars->SetNumberOfComponents(1);
    cellScalars->Allocate(sf.size());
    cellScalars->SetName(sf.name().c_str());

    forAll(superCells_, sci)
    {
        cellScalars->InsertComponent(sci, 0, sf[superCells_[sci]]);
    }

    vtkMesh->GetCellData()->AddArray(cellScalars);
    if (!vtkMesh->GetCellData()->GetScalars())
    {
        vtkMesh->GetCellData()->SetScalars(cellScalars);
    }

    cellScalars->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
