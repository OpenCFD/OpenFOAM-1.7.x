/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#include "backwardsCompatibilityWallFunctions.H"

#include "calculatedFvPatchField.H"
#include "alphatWallFunctionFvPatchScalarField.H"
#include "mutWallFunctionFvPatchScalarField.H"
#include "epsilonWallFunctionFvPatchScalarField.H"
#include "kqRWallFunctionFvPatchField.H"
#include "omegaWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<volScalarField> autoCreateAlphat
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    IOobject alphatHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (alphatHeader.headerOk())
    {
        return tmp<volScalarField>(new volScalarField(alphatHeader, mesh));
    }
    else
    {
        Info<< "--> Creating " << fieldName
            << " to employ run-time selectable wall functions" << endl;

        const fvBoundaryMesh& bm = mesh.boundary();

        wordList alphatBoundaryTypes(bm.size());

        forAll(bm, patchI)
        {
            if (isA<wallFvPatch>(bm[patchI]))
            {
                alphatBoundaryTypes[patchI] =
                    RASModels::alphatWallFunctionFvPatchScalarField::typeName;
            }
            else
            {
                alphatBoundaryTypes[patchI] =
                    calculatedFvPatchField<scalar>::typeName;
            }
        }

        tmp<volScalarField> alphat
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimDensity*dimArea/dimTime, 0.0),
                alphatBoundaryTypes
            )
        );

        Info<< "    Writing new " << fieldName << endl;
        alphat().write();

        return alphat;
    }
}


tmp<volScalarField> autoCreateMut
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    IOobject mutHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (mutHeader.headerOk())
    {
        return tmp<volScalarField>(new volScalarField(mutHeader, mesh));
    }
    else
    {
        Info<< "--> Creating " << fieldName
            << " to employ run-time selectable wall functions" << endl;

        const fvBoundaryMesh& bm = mesh.boundary();

        wordList mutBoundaryTypes(bm.size());

        forAll(bm, patchI)
        {
            if (isA<wallFvPatch>(bm[patchI]))
            {
                mutBoundaryTypes[patchI] =
                    RASModels::mutWallFunctionFvPatchScalarField::typeName;
            }
            else
            {
                mutBoundaryTypes[patchI] =
                    calculatedFvPatchField<scalar>::typeName;
            }
        }

        tmp<volScalarField> mut
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimensionedScalar("zero", dimDensity*dimArea/dimTime, 0.0),
                mutBoundaryTypes
            )
        );

        Info<< "    Writing new " << fieldName << endl;
        mut().write();

        return mut;
    }
}


tmp<volScalarField> autoCreateEpsilon
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::epsilonWallFunctionFvPatchScalarField
        >
        (
            fieldName,
            mesh
        );
}


tmp<volScalarField> autoCreateOmega
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::omegaWallFunctionFvPatchScalarField
        >
        (
            fieldName,
            mesh
        );
}


tmp<volScalarField> autoCreateK
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::kqRWallFunctionFvPatchField<scalar>
        >
        (
            fieldName,
            mesh
        );
}


tmp<volScalarField> autoCreateQ
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            scalar,
            RASModels::kqRWallFunctionFvPatchField<scalar>
        >
        (
            fieldName,
            mesh
        );
}


tmp<volSymmTensorField> autoCreateR
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return
        autoCreateWallFunctionField
        <
            symmTensor,
            RASModels::kqRWallFunctionFvPatchField<symmTensor>
        >
        (
            fieldName,
            mesh
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //

