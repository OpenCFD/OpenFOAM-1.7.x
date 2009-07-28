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
    Dimensioned sphericalTensor obtained from generic dimensioned type.

\*---------------------------------------------------------------------------*/

#include "dimensionedSphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
dimensionedSphericalTensor dimensionedSphericalTensor::T() const
{
    return dimensionedSphericalTensor
    (
        name()+".T()",
        dimensions(),
        value().T()
    );
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

dimensionedScalar tr(const dimensionedSphericalTensor& dt)
{
    return dimensionedScalar
    (
        "tr("+dt.name()+')',
        dt.dimensions(),
        tr(dt.value())
    );
}


dimensionedScalar det(const dimensionedSphericalTensor& dt)
{
    return dimensionedScalar
    (
        "det("+dt.name()+')',
        pow(dt.dimensions(), sphericalTensor::dim),
        det(dt.value())
    );
}


dimensionedSphericalTensor inv(const dimensionedSphericalTensor& dt)
{
    return dimensionedSphericalTensor
    (
        "inv("+dt.name()+')',
        dimless/dt.dimensions(),
        inv(dt.value())
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
