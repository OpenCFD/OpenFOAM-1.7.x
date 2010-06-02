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

#include "dynMixedSmagorinsky.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynMixedSmagorinsky, 0);
addToRunTimeSelectionTable(LESModel, dynMixedSmagorinsky, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynMixedSmagorinsky::dynMixedSmagorinsky
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    LESModel(typeName, U, phi, transport),
    scaleSimilarity(U, phi, transport),
    dynSmagorinsky(U, phi, transport)
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> dynMixedSmagorinsky::k() const
{
    return
    (
        scaleSimilarity::k()
      + dynSmagorinsky::k()
    );
}


tmp<volScalarField> dynMixedSmagorinsky::epsilon() const
{
    return
    (
        scaleSimilarity::epsilon()
      + dynSmagorinsky::epsilon()
    );
}


tmp<volSymmTensorField> dynMixedSmagorinsky::B() const
{
    return
    (
        scaleSimilarity::B()
      + dynSmagorinsky::B()
    );
}


tmp<volSymmTensorField> dynMixedSmagorinsky::devBeff() const
{
    return
    (
        scaleSimilarity::devBeff()
      + dynSmagorinsky::devBeff()
    );
}


tmp<fvVectorMatrix> dynMixedSmagorinsky::divDevBeff(volVectorField& U) const
{
    return
    (
        scaleSimilarity::divDevBeff(U)
      + dynSmagorinsky::divDevBeff(U)
    );
}


void dynMixedSmagorinsky::correct(const tmp<volTensorField>& gradU)
{
    scaleSimilarity::correct(gradU);
    dynSmagorinsky::correct(gradU());
}


bool dynMixedSmagorinsky::read()
{
    if (LESModel::read())
    {
        scaleSimilarity::read();
        dynSmagorinsky::read();

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // namespace incompressible
} // End namespace Foam

// ************************************************************************* //
