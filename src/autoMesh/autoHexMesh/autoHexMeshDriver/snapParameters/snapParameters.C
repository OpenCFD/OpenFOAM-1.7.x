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

#include "snapParameters.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::snapParameters::snapParameters(const dictionary& dict, const label dummy)
:
    nSmoothPatch_(readLabel(dict.lookup("nSmoothPatch"))),
    snapTol_(readScalar(dict.lookup("snapTol"))),
    nSmoothDispl_(readLabel(dict.lookup("nSmoothDispl"))),
    nSnap_(readLabel(dict.lookup("nSnap")))
{}


// Construct from dictionary
Foam::snapParameters::snapParameters(const dictionary& dict)
:
    nSmoothPatch_(readLabel(dict.lookup("nSmoothPatch"))),
    snapTol_(readScalar(dict.lookup("tolerance"))),
    nSmoothDispl_(readLabel(dict.lookup("nSolveIter"))),
    nSnap_(readLabel(dict.lookup("nRelaxIter")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
