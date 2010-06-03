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

#include "laminarFlameSpeed.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::laminarFlameSpeed> Foam::laminarFlameSpeed::New
(
    const hhuCombustionThermo& ct
)
{
    IOdictionary laminarFlameSpeedDict
    (
        IOobject
        (
            "combustionProperties",
            ct.T().time().constant(),
            ct.T().db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word laminarFlameSpeedType
    (
        laminarFlameSpeedDict.lookup("laminarFlameSpeedCorrelation")
    );

    Info<< "Selecting laminar flame speed correlation "
        << laminarFlameSpeedType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(laminarFlameSpeedType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "laminarFlameSpeed::New(const hhuCombustionThermo&)",
            laminarFlameSpeedDict
        )   << "Unknown laminarFlameSpeed type "
            << laminarFlameSpeedType << endl << endl
            << "Valid laminarFlameSpeed types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<laminarFlameSpeed>(cstrIter()(laminarFlameSpeedDict, ct));
}


// ************************************************************************* //
