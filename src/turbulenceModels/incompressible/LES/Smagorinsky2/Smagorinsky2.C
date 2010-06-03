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

#include "Smagorinsky2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Smagorinsky2, 0);
addToRunTimeSelectionTable(LESModel, Smagorinsky2, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Smagorinsky2::Smagorinsky2
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
:
    LESModel(typeName, U, phi, transport),
    Smagorinsky(U, phi, transport),

    cD2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cD2",
            coeffDict_,
            0.02
        )
    )
{
    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Evaluate B (from the definition of an eddy viscosity model) and
// return it.

tmp<volSymmTensorField> Smagorinsky2::B() const
{
    volSymmTensorField D = dev(symm(fvc::grad(U())));

    return (((2.0/3.0)*I)*k() - 2.0*nuSgs_*D - (2.0*cD2_)*delta()*(D&D));
}


tmp<fvVectorMatrix> Smagorinsky2::divDevBeff
(
    volVectorField& U
) const
{
    volTensorField gradU = fvc::grad(U);

    volSymmTensorField aniNuEff
    (
        "aniNuEff",
        I*nuEff() + cD2_*delta()*symm(gradU)
    );

    return
    (
      - fvm::laplacian(aniNuEff, U) - fvc::div(nuEff()*dev(fvc::grad(U)().T()))
    );
}


bool Smagorinsky2::read()
{
    if (Smagorinsky::read())
    {
        cD2.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
