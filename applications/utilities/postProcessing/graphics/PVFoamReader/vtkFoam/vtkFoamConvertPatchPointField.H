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

#ifndef vtkFoamConvertPatchPointField_H
#define vtkFoamConvertPatchPointField_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::vtkFoam::convertPatchPointField
(
    const word& name,
    const Field<Type>& pptf,
    const label regioni
)
{
    vtkUnstructuredGrid *vtkMesh =
        vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(regioni));

    vtkFloatArray *pointTypes = vtkFloatArray::New();
    pointTypes->SetNumberOfTuples(pptf.size());
    pointTypes->SetNumberOfComponents(Type::nComponents);
    pointTypes->Allocate(Type::nComponents*pptf.size());
    pointTypes->SetName(name.c_str());

    float vec[Type::nComponents];

    forAll(pptf, i)
    {
        for (direction d=0; d<Type::nComponents; d++)
        {
            vec[d] = pptf[i][d];
        }

        pointTypes->InsertTuple(i, vec);
    }

    vtkMesh->GetPointData()->AddArray(pointTypes);
    pointTypes->Delete();
}


template<>
void Foam::vtkFoam::convertPatchPointField
(
    const word& name,
    const Field<scalar>& ppsf,
    const label regioni
)
{
    vtkUnstructuredGrid *vtkMesh =
        vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(regioni));

    vtkFloatArray *pointScalars = vtkFloatArray::New();
    pointScalars->SetNumberOfTuples(ppsf.size());
    pointScalars->SetNumberOfComponents(1);
    pointScalars->Allocate(ppsf.size());
    pointScalars->SetName(name.c_str());

    for (int i=0; i<ppsf.size(); i++)
    {
        pointScalars->InsertComponent(i, 0, ppsf[i]);
    }

    vtkMesh->GetPointData()->AddArray(pointScalars);
    if (!vtkMesh->GetPointData()->GetScalars())
    {
        vtkMesh->GetPointData()->SetScalars(pointScalars);
    }

    pointScalars->Delete();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
