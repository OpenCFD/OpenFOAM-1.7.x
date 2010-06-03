/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
inline Foam::PtrList<Foam::scalarField>&
Foam::ODEChemistryModel<CompType, ThermoType>::RR()
{
    return RR_;
}


template<class CompType, class ThermoType>
inline const Foam::PtrList<Foam::Reaction<ThermoType> >&
Foam::ODEChemistryModel<CompType, ThermoType>::reactions() const
{
    return reactions_;
}


template<class CompType, class ThermoType>
inline const Foam::PtrList<ThermoType>&
Foam::ODEChemistryModel<CompType, ThermoType>::specieThermo() const
{
    return specieThermo_;
}


template<class CompType, class ThermoType>
inline Foam::label
Foam::ODEChemistryModel<CompType, ThermoType>::nSpecie() const
{
    return nSpecie_;
}


template<class CompType, class ThermoType>
inline Foam::label
Foam::ODEChemistryModel<CompType, ThermoType>::nReaction() const
{
    return nReaction_;
}


template<class CompType, class ThermoType>
inline const Foam::chemistrySolver<CompType, ThermoType>&
Foam::ODEChemistryModel<CompType, ThermoType>::solver() const
{
    return solver_;
}


template<class CompType, class ThermoType>
inline Foam::tmp<Foam::volScalarField>
Foam::ODEChemistryModel<CompType, ThermoType>::RR
(
    const label i
) const
{
    tmp<volScalarField> tRR
    (
        new volScalarField
        (
            IOobject
            (
                "RR(" + this->Y_[i].name() + ')',
                this->time().timeName(),
                this->mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->mesh(),
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (this->chemistry_)
    {
        tRR().internalField() = RR_[i];
        tRR().correctBoundaryConditions();
    }
    return tRR;
}


// ************************************************************************* //
