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

InClass
    Foam::vtkFoam

\*---------------------------------------------------------------------------*/

#ifndef vtkFoamConvertPatchFaceField_H
#define vtkFoamConvertPatchFaceField_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkFoam::convertPatchFaceField
(
    const word& name,
    const Field<Type>& ptf,
    const label regioni
)
{
    vtkUnstructuredGrid *vtkMesh =
        vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(regioni));

    vtkFloatArray *cellTypes = vtkFloatArray::New();
    cellTypes->SetNumberOfTuples(ptf.size());
    cellTypes->SetNumberOfComponents(Type::nComponents);
    cellTypes->Allocate(Type::nComponents*ptf.size());
    cellTypes->SetName(name.c_str());

    float vec[Type::nComponents];

    forAll(ptf, i)
    {
        const Type& t = ptf[i];
        for (direction d=0; d<Type::nComponents; d++)
        {
            vec[d] = t[d];
        }

        cellTypes->InsertTuple(i, vec);
    }

    vtkMesh->GetCellData()->AddArray(cellTypes);
    cellTypes->Delete();
}


template<>
void Foam::vtkFoam::convertPatchFaceField
(
    const word& name,
    const Field<scalar>& psf,
    const label regioni
)
{
    vtkUnstructuredGrid *vtkMesh =
        vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(regioni));

    vtkFloatArray *cellScalars = vtkFloatArray::New();
    cellScalars->SetNumberOfTuples(psf.size());
    cellScalars->SetNumberOfComponents(1);
    cellScalars->Allocate(psf.size());
    cellScalars->SetName(name.c_str());

    forAll(psf, i)
    {
        cellScalars->InsertComponent(i, 0, psf[i]);
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
