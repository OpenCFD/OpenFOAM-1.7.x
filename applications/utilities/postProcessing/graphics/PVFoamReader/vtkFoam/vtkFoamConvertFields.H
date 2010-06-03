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

#ifndef vtkFoamConvertFields_H
#define vtkFoamConvertFields_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "vtkFoamConvertVolField.H"
#include "vtkFoamConvertPointField.H"

#include "vtkFoamConvertPatchFaceField.H"
#include "vtkFoamConvertPatchPointField.H"

#include "emptyFvPatchField.H"

template<class Type>
void Foam::vtkFoam::convertVolFields
(
    const fvMesh& mesh,
    const volPointInterpolation& pInterp,
    const IOobjectList& objects,
    vtkDataArraySelection *fieldSelection
)
{
    IOobjectList fieldObjects
    (
        objects.lookupClass
        (
            GeometricField<Type, fvPatchField, volMesh>::typeName
        )
    );

    label nSelectedFields = fieldSelection->GetNumberOfArrays();

    for (label i=0; i<nSelectedFields; i++)
    {
        if(fieldSelection->GetArraySetting(i))
        {
            word fieldName = fieldSelection->GetArrayName(i);

            if (fieldObjects.found(fieldName))
            {
                GeometricField<Type, fvPatchField, volMesh> tf
                (
                    IOobject
                    (
                        fieldName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                );

                tmp<GeometricField<Type, pointPatchField, pointMesh> > tptf
                (
                    pInterp.interpolate(tf)
                );

                if (selectedRegions_[0])
                {
                    convertVolField(tf);
                    convertPointField(tptf(), tf);
                }

                label regioni = 0;
                forAll (mesh.boundaryMesh(), patchi)
                {
                    if (mesh.boundaryMesh()[patchi].size())
                    {
                        regioni++;

                        if (selectedRegions_[regioni])
                        {
                            const fvPatchField<Type>& ptf
                            (
                                tf.boundaryField()[patchi]
                            );

                            if (!isType<emptyFvPatchField<Type> >(ptf))
                            {
                                convertPatchFaceField
                                (
                                    tf.name(),
                                    ptf,
                                    regioni
                                );

                                convertPatchPointField
                                (
                                    tptf().name(),
                                    tptf().boundaryField()[patchi]
                                        .patchInternalField()(),
                                    regioni
                                );
                            }
                            else
                            {
                                fvPatch p
                                (
                                    ptf.patch().patch(),
                                    tf.mesh().boundary()
                                );

                                convertPatchFaceField
                                (
                                    tf.name(),
                                    fvPatchField<Type>(p, tf)
                                        .patchInternalField()(),
                                    regioni
                                );

                                convertPatchPointField
                                (
                                    tptf().name(),
                                    tptf().boundaryField()[patchi]
                                       .patchInternalField()(),
                                    regioni
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::vtkFoam::convertPointFields
(
    const fvMesh& mesh,
    const IOobjectList& objects,
    vtkDataArraySelection *fieldSelection
)
{
    IOobjectList fieldObjects
    (
        objects.lookupClass
        (
            GeometricField<Type, pointPatchField, pointMesh>::typeName
        )
    );

    label nSelectedFields = fieldSelection->GetNumberOfArrays();

    for (label i=0; i<nSelectedFields; i++)
    {
        if(fieldSelection->GetArraySetting(i))
        {
            word fieldName = fieldSelection->GetArrayName(i);

            if (fieldObjects.found(fieldName))
            {
                pointMesh pMesh(mesh);

                GeometricField<Type, pointPatchField, pointMesh> ptf
                (
                    IOobject
                    (
                        fieldName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    pMesh
                );

                if (selectedRegions_[0])
                {
                    convertPointField
                    (
                        ptf,
                        GeometricField<Type, fvPatchField, volMesh>::null()
                    );
                }

                label regioni = 0;
                forAll (mesh.boundaryMesh(), patchi)
                {
                    if (mesh.boundaryMesh()[patchi].size())
                    {
                        regioni++;

                        if (selectedRegions_[regioni])
                        {
                            convertPatchPointField
                            (
                                ptf.name(),
                                ptf.boundaryField()[patchi]
                               .patchInternalField()(),
                                regioni
                            );
                        }
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
