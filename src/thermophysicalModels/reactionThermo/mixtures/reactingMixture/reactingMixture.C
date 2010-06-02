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

#include "reactingMixture.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::reactingMixture<ThermoType>::reactingMixture
(
    const dictionary& thermoDict,
    const fvMesh& mesh
)
:
    autoPtr<chemistryReader<ThermoType> >
    (
        chemistryReader<ThermoType>::New(thermoDict)
    ),
    multiComponentMixture<ThermoType>
    (
        thermoDict,
        autoPtr<chemistryReader<ThermoType> >::operator()().species(),
        autoPtr<chemistryReader<ThermoType> >::operator()().speciesThermo(),
        mesh
    ),
    PtrList<Reaction<ThermoType> >
    (
        autoPtr<chemistryReader<ThermoType> >::operator()().reactions(),
        this->species_
    )
{
    autoPtr<chemistryReader<ThermoType> >::clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::reactingMixture<ThermoType>::read(const dictionary& thermoDict)
{}


// ************************************************************************* //
