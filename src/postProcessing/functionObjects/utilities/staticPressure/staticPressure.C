/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "staticPressure.H"
#include "volFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(staticPressure, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::staticPressure::isKinematicPressure()
{
    const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

    return p.dimensions() == sqr(dimLength)/sqr(dimTime);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::staticPressure::staticPressure
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    rho_(readScalar(dict.lookup("rho")))
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "staticPressure::staticPressure"
            "(const objectRegistry&, const dictionary&)"
        )   << "No fvMesh available, deactivating." << nl
            << endl;
    }
    else
    {
        // Check if the pressure is kinematic pressure, otherwise deactivate
        if (!isKinematicPressure())
        {
            active_ = false;
            WarningIn
            (
                "staticPressure::staticPressure"
                "(const objectRegistry&, const dictionary&)"
            )   << "Pressure is not kinematic pressure, deactivating." << nl
                << endl;
        }
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::staticPressure::~staticPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::staticPressure::read(const dictionary& dict)
{
    if (active_)
    {
        dict.readIfPresent("p", pName_);
        dict.lookup("rho") >> rho_;
    }
}


void Foam::staticPressure::execute()
{
    // Do nothing - only valid on write
}


void Foam::staticPressure::end()
{
    // Do nothing - only valid on write
}


void Foam::staticPressure::write()
{
    if (active_)
    {
        const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

        volScalarField pStatic
        (
            IOobject
            (
                "pStatic",
                obr_.time().timeName(),
                obr_,
                IOobject::NO_READ
            ),
            dimensionedScalar("rho", dimDensity, rho_)*p
        );

        pStatic.write();
    }
}


// ************************************************************************* //
