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

\*---------------------------------------------------------------------------*/

#include "fvcAverage.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMesh.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
average
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    const fvMesh& mesh = ssf.mesh();

    tmp<GeometricField<Type, fvPatchField, volMesh> > taverage
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "average("+ssf.name()+')',
                ssf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            ssf.dimensions()
        )
    );

    GeometricField<Type, fvPatchField, volMesh>& av = taverage();

    av.internalField() = 
    (
        surfaceSum(mesh.magSf()*ssf)/surfaceSum(mesh.magSf())
    )().internalField();

    forAll(av.boundaryField(), patchi)
    {
        av.boundaryField()[patchi] = ssf.boundaryField()[patchi];
    }

    av.correctBoundaryConditions();

    return taverage;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
average
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >& tssf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > taverage
    (
        fvc::average(tssf())
    );
    tssf.clear();
    return taverage;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
average
(
    const GeometricField<Type, fvPatchField, volMesh>& vtf
)
{
    return fvc::average(linearInterpolate(vtf));
}


template<class Type> 
tmp<GeometricField<Type, fvPatchField, volMesh> >
average
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvtf
)
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > taverage
    (
        fvc::average(tvtf())
    );
    tvtf.clear();
    return taverage;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
