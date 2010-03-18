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

\*----------------------------------------------------------------------------*/

#include "porousZone.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

// adjust negative resistance values to be multiplier of max value
void Foam::porousZone::adjustNegativeResistance(dimensionedVector& resist)
{
    scalar maxCmpt = max(0, cmptMax(resist.value()));

    if (maxCmpt < 0)
    {
        FatalErrorIn
        (
            "Foam::porousZone::porousZone::adjustNegativeResistance"
            "(dimensionedVector&)"
        )   << "negative resistances! " << resist
            << exit(FatalError);
    }
    else
    {
        vector& val = resist.value();
        for (label cmpt=0; cmpt < vector::nComponents; ++cmpt)
        {
            if (val[cmpt] < 0)
            {
                val[cmpt] *= -maxCmpt;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousZone::porousZone
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    cellZoneID_(mesh_.cellZones().findZoneID(name)),
    coordSys_(dict, mesh),
    porosity_(1),
    C0_(0),
    C1_(0),
    D_("D", dimensionSet(0, -2, 0, 0, 0), tensor::zero),
    F_("F", dimensionSet(0, -1, 0, 0, 0), tensor::zero)
{
    Info<< "Creating porous zone: " << name_ << endl;

    bool foundZone = (cellZoneID_ != -1);
    reduce(foundZone, orOp<bool>());

    if (!foundZone && Pstream::master())
    {
        FatalErrorIn
        (
            "Foam::porousZone::porousZone"
            "(const fvMesh&, const word&, const dictionary&)"
        )   << "cannot find porous cellZone " << name_
            << exit(FatalError);
    }


    // porosity
    if (dict_.readIfPresent("porosity", porosity_))
    {
        if (porosity_ <= 0.0 || porosity_ > 1.0)
        {
            FatalIOErrorIn
            (
                "Foam::porousZone::porousZone"
                "(const fvMesh&, const word&, const dictionary&)",
                dict_
            )
                << "out-of-range porosity value " << porosity_
                << exit(FatalIOError);
        }
    }

    // powerLaw coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("powerLaw"))
    {
        dictPtr->readIfPresent("C0", C0_);
        dictPtr->readIfPresent("C1", C1_);
    }

    // Darcy-Forchheimer coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("Darcy"))
    {
        // local-to-global transformation tensor
        const tensor& E = coordSys_.R();

        dimensionedVector d(vector::zero);
        if (dictPtr->readIfPresent("d", d))
        {
            if (D_.dimensions() != d.dimensions())
            {
                FatalIOErrorIn
                (
                    "Foam::porousZone::porousZone"
                    "(const fvMesh&, const word&, const dictionary&)",
                    dict_
                )   << "incorrect dimensions for d: " << d.dimensions()
                    << " should be " << D_.dimensions()
                    << exit(FatalIOError);
            }

            adjustNegativeResistance(d);

            D_.value().xx() = d.value().x();
            D_.value().yy() = d.value().y();
            D_.value().zz() = d.value().z();
            D_.value() = (E & D_ & E.T()).value();
        }

        dimensionedVector f(vector::zero);
        if (dictPtr->readIfPresent("f", f))
        {
            if (F_.dimensions() != f.dimensions())
            {
                FatalIOErrorIn
                (
                    "Foam::porousZone::porousZone"
                    "(const fvMesh&, const word&, const dictionary&)",
                    dict_
                )   << "incorrect dimensions for f: " << f.dimensions()
                    << " should be " << F_.dimensions()
                    << exit(FatalIOError);
            }

            adjustNegativeResistance(f);

            // leading 0.5 is from 1/2 * rho
            F_.value().xx() = 0.5*f.value().x();
            F_.value().yy() = 0.5*f.value().y();
            F_.value().zz() = 0.5*f.value().z();
            F_.value() = (E & F_ & E.T()).value();
        }
    }

    // provide some feedback for the user
    // writeDict(Info, false);

    // it is an error not to define anything
    if
    (
        C0_ <= VSMALL
     && magSqr(D_.value()) <= VSMALL
     && magSqr(F_.value()) <= VSMALL
    )
    {
        FatalIOErrorIn
        (
            "Foam::porousZone::porousZone"
            "(const fvMesh&, const word&, const dictionary&)",
            dict_
        )   << "neither powerLaw (C0/C1) "
               "nor Darcy-Forchheimer law (d/f) specified"
            << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousZone::addResistance(fvVectorMatrix& UEqn) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    bool compressible = false;
    if (UEqn.dimensions() == dimensionSet(1, 1, -2, 0, 0))
    {
        compressible = true;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();

    if (C0_ > VSMALL)
    {
        if (compressible)
        {
            addPowerLawResistance
            (
                Udiag,
                cells,
                V,
                mesh_.lookupObject<volScalarField>("rho"),
                U
            );
        }
        else
        {
            addPowerLawResistance
            (
                Udiag,
                cells,
                V,
                geometricOneField(),
                U
            );
        }
    }

    const tensor& D = D_.value();
    const tensor& F = F_.value();

    if (magSqr(D) > VSMALL || magSqr(F) > VSMALL)
    {
        if (compressible)
        {
            addViscousInertialResistance
            (
                Udiag,
                Usource,
                cells,
                V,
                mesh_.lookupObject<volScalarField>("rho"),
                mesh_.lookupObject<volScalarField>("mu"),
                U
            );
        }
        else
        {
            addViscousInertialResistance
            (
                Udiag,
                Usource,
                cells,
                V,
                geometricOneField(),
                mesh_.lookupObject<volScalarField>("nu"),
                U
            );
        }
    }
}


void Foam::porousZone::addResistance
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU,
    bool correctAUprocBC
) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    bool compressible = false;
    if (UEqn.dimensions() == dimensionSet(1, 1, -2, 0, 0))
    {
        compressible = true;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const vectorField& U = UEqn.psi();

    if (C0_ > VSMALL)
    {
        if (compressible)
        {
            addPowerLawResistance
            (
                AU,
                cells,
                mesh_.lookupObject<volScalarField>("rho"),
                U
            );
        }
        else
        {
            addPowerLawResistance
            (
                AU,
                cells,
                geometricOneField(),
                U
            );
        }
    }

    const tensor& D = D_.value();
    const tensor& F = F_.value();

    if (magSqr(D) > VSMALL || magSqr(F) > VSMALL)
    {
        if (compressible)
        {
            addViscousInertialResistance
            (
                AU,
                cells,
                mesh_.lookupObject<volScalarField>("rho"),
                mesh_.lookupObject<volScalarField>("mu"),
                U
            );
        }
        else
        {
            addViscousInertialResistance
            (
                AU,
                cells,
                geometricOneField(),
                mesh_.lookupObject<volScalarField>("nu"),
                U
            );
        }
    }

    if (correctAUprocBC)
    {
        // Correct the boundary conditions of the tensorial diagonal to ensure
        // processor boundaries are correctly handled when AU^-1 is interpolated
        // for the pressure equation.
        AU.correctBoundaryConditions();
    }
}


void Foam::porousZone::writeDict(Ostream& os, bool subDict) const
{
    if (subDict)
    {
        os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
        os.writeKeyword("name") << zoneName() << token::END_STATEMENT << nl;
    }
    else
    {
        os  << indent << zoneName() << nl
            << indent << token::BEGIN_BLOCK << incrIndent << nl;
    }

    if (dict_.found("note"))
    {
        os.writeKeyword("note") << string(dict_.lookup("note"))
            << token::END_STATEMENT << nl;
    }

    coordSys_.writeDict(os, true);

    if (dict_.found("porosity"))
    {
        os.writeKeyword("porosity") << porosity() << token::END_STATEMENT << nl;
    }

    // powerLaw coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("powerLaw"))
    {
        os << indent << "powerLaw";
        dictPtr->write(os);
    }

    // Darcy-Forchheimer coefficients
    if (const dictionary* dictPtr = dict_.subDictPtr("Darcy"))
    {
        os << indent << "Darcy";
        dictPtr->write(os);
    }

    os << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const porousZone& pZone)
{
    pZone.writeDict(os);
    return os;
}

// ************************************************************************* //
