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

InClass
    Foam::vtkFoam

\*---------------------------------------------------------------------------*/

#ifndef vtkFoamConvertPointField_H
#define vtkFoamConvertPointField_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interpolatePointToCell.H"

template<class Type>
void Foam::vtkFoam::convertPointField
(
    const GeometricField<Type, pointPatchField, pointMesh>& ptf,
    const GeometricField<Type, fvPatchField, volMesh>& tf
)
{
    vtkUnstructuredGrid *vtkMesh =
        vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(0));

    vtkFloatArray *pointTypes = vtkFloatArray::New();
    pointTypes->SetNumberOfTuples(ptf.size() + addPointCellLabels_.size());
    pointTypes->SetNumberOfComponents(Type::nComponents);
    pointTypes->Allocate(Type::nComponents*ptf.size());
    pointTypes->SetName(ptf.name().c_str());

    float vec[Type::nComponents];

    forAll(ptf, i)
    {
        for (direction d=0; d<Type::nComponents; d++)
        {
            vec[d] = ptf[i][d];
        }

        pointTypes->InsertTuple(i, vec);
    }

    label i = ptf.size();

    if (&tf != &GeometricField<Type, fvPatchField, volMesh>::null())
    {
        forAll(addPointCellLabels_, api)
        {
            Type t = tf[addPointCellLabels_[api]];

            for (direction d=0; d<Type::nComponents; d++)
            {
                vec[d] = t[d];
            }

            pointTypes->InsertTuple(i++, vec);
        }
    }
    else
    {
        forAll(addPointCellLabels_, api)
        {
            Type t = interpolatePointToCell(ptf, addPointCellLabels_[api]);

            for (direction d=0; d<Type::nComponents; d++)
            {
                vec[d] = t[d];
            }

            pointTypes->InsertTuple(i++, vec);
        }
    }

    vtkMesh->GetPointData()->AddArray(pointTypes);
    pointTypes->Delete();
}


template<>
void Foam::vtkFoam::convertPointField
(
    const GeometricField<scalar, pointPatchField, pointMesh>& psf,
    const GeometricField<scalar, fvPatchField, volMesh>& sf
)
{
    vtkUnstructuredGrid *vtkMesh =
        vtkUnstructuredGrid::SafeDownCast(reader_->GetOutput(0));

    vtkFloatArray *pointScalars = vtkFloatArray::New();
    pointScalars->SetNumberOfTuples(psf.size() + addPointCellLabels_.size());
    pointScalars->SetNumberOfComponents(1);
    pointScalars->Allocate(psf.size());
    pointScalars->SetName(psf.name().c_str());

    for (int i=0; i<psf.size(); i++)
    {
        pointScalars->InsertComponent(i, 0, psf[i]);
    }

    label i = psf.size();

    if (&sf != &GeometricField<scalar, fvPatchField, volMesh>::null())
    {
        forAll(addPointCellLabels_, api)
        {
            pointScalars->InsertComponent
            (
                i++,
                0,
                sf[addPointCellLabels_[api]]
            );
        }
    }
    else
    {
        forAll(addPointCellLabels_, api)
        {
            pointScalars->InsertComponent
            (
                i++,
                0,
                interpolatePointToCell(psf, addPointCellLabels_[api])
            );
        }
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
