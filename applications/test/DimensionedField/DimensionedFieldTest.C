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

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "DimensionedFields.H"
#include "DimensionedSphericalTensorField.H"
#include "vector.H"
#include "tensor.H"
#include "GeoMesh.H"

using namespace Foam;

namespace Foam
{

class vMesh
{

public:

    vMesh()
    {}

    label size() const
    {
        return 10;
    }
};

};

template<>
const word Foam::DimensionedField<scalar, GeoMesh<vMesh> >::typeName
(
    "dimenionedScalarField"
);

template<>
const word Foam::DimensionedField<vector, GeoMesh<vMesh> >::typeName
(
    "dimenionedVectorField"
);

template<>
const word Foam::DimensionedField<tensor, GeoMesh<vMesh> >::typeName
(
    "dimenionedTensorField"
);

template<>
const word Foam::DimensionedField<sphericalTensor, GeoMesh<vMesh> >::typeName
(
    "dimenionedSphericalTensorField"
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"

    vMesh vm;

    DimensionedField<scalar, GeoMesh<vMesh> > dsf
    (
        IOobject
        (
            "dsf",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        vm
    );

    Info<< dsf << endl;
    dsf += dsf;
    dsf -= dimensionedScalar("5", dsf.dimensions(), 5.0);
    Info<< dsf << endl;

    Info<< sqr(dsf + dsf) - sqr(dsf + dsf) << endl;

    DimensionedField<vector, GeoMesh<vMesh> > dvf
    (
        IOobject
        (
            "dvf",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        vm
    );

    Info<< (dvf ^ (dvf ^ dvf)) << endl;

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
