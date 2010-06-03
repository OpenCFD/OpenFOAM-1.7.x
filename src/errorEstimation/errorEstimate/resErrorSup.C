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

\*---------------------------------------------------------------------------*/

#include "resErrorSup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace resError
{

template<class Type>
tmp<errorEstimate<Type> >
Sp
(
    const volScalarField& sp,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<errorEstimate<Type> >
    (
        new errorEstimate<Type>
        (
            vf,
            sp.dimensions()*vf.dimensions(),
            sp.internalField()*vf.internalField(),
            scalarField(vf.internalField().size(), 0)
        )
    );
}

template<class Type>
tmp<errorEstimate<Type> >
Sp
(
    const tmp<volScalarField>& tsp,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<errorEstimate<Type> > tee = resError::Sp(tsp(), vf);
    tsp.clear();
    return tee;
}


template<class Type>
tmp<errorEstimate<Type> >
Sp
(
    const dimensionedScalar& sp,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return tmp<errorEstimate<Type> >
    (
        new errorEstimate<Type>
        (
            vf,
            sp.dimensions()*vf.dimensions(),
            sp.value()*vf.internalField(),
            scalarField(vf.internalField().size(), 0)
        )
    );
}


template<class Type>
tmp<errorEstimate<Type> >
SuSp
(
    const volScalarField& sp,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return Sp(sp, vf);
}

template<class Type>
tmp<errorEstimate<Type> >
SuSp
(
    const tmp<volScalarField>& tsp,
    GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<errorEstimate<Type> > tee = resError::SuSp(tsp(), vf);
    tsp.clear();
    return tee;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace resError

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

