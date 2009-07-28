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

#include "solidBodyMotionFunction.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidBodyMotionFunction> Foam::solidBodyMotionFunction::New
(
    const dictionary& SBMFCoeffs,
    const Time& runTime
)
{
    word solidBodyMotionFunctionTypeName = 
        SBMFCoeffs.lookup("solidBodyMotionFunction");

    Info<< "Selecting solid-body motion function "
        << solidBodyMotionFunctionTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solidBodyMotionFunctionTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "solidBodyMotionFunction::New"
            "("
            "    const dictionary& SBMFCoeffs,"
            "    const Time& runTime"
            ")"
        )   << "Unknown solidBodyMotionFunction type "
            << solidBodyMotionFunctionTypeName << endl << endl
            << "Valid  solidBodyMotionFunctions are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<solidBodyMotionFunction>(cstrIter()(SBMFCoeffs, runTime));
}


// ************************************************************************* //
