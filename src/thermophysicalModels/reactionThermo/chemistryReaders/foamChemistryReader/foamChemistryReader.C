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

#include "foamChemistryReader.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::foamChemistryReader<ThermoType>::foamChemistryReader
(
    const fileName& reactionsFileName,
    const fileName& thermoFileName
)
:
    chemistryReader<ThermoType>(),
    speciesThermo_(IFstream(thermoFileName)()),
    speciesTable_(dictionary(IFstream(reactionsFileName)()).lookup("species")),
    reactions_
    (
        dictionary(IFstream(reactionsFileName)()).lookup("reactions"),
        Reaction<ThermoType>::iNew(speciesTable_, speciesThermo_)
    )
{}


template<class ThermoType>
Foam::foamChemistryReader<ThermoType>::foamChemistryReader
(
    const dictionary& thermoDict
)
:
    chemistryReader<ThermoType>(),
    speciesThermo_
    (
        IFstream
        (
            fileName(thermoDict.lookup("foamChemistryThermoFile")).expand()
        )()
    ),
    speciesTable_
    (
        dictionary
        (
            IFstream
            (
                fileName(thermoDict.lookup("foamChemistryFile")).expand()
            )()
        ).lookup("species")
    ),
    reactions_
    (
        dictionary
        (
            IFstream
            (
                fileName(thermoDict.lookup("foamChemistryFile")).expand()
            )()
        ).lookup("reactions"),
        typename Reaction<ThermoType>::iNew(speciesTable_, speciesThermo_)
    )
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
