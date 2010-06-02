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

Description
    Miscellaneous collection of functions and template related to Ensight data

SourceFiles
    ensightOutputFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef ensightOutputFunctions_H
#define ensightOutputFunctions_H

#include "ensightFile.H"
#include "Cloud.H"
#include "polyMesh.H"
#include "IOobject.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void ensightCaseEntry
(
    OFstream& caseFile,
    const string& ensightType,
    const word& fieldName,
    const fileName& dataMask,
    const fileName& local=fileName::null,
    const label cloudNo=-1,
    const label timeSet=1
);


void ensightParticlePositions
(
    const polyMesh& mesh,
    const fileName& dataDir,
    const fileName& subDir,
    const word& cloudName,
    IOstream::streamFormat format
);


//- write lagrangian parcels
template<class Type>
void ensightLagrangianField
(
    const IOobject& fieldObject,
    const fileName& dataDir,
    const fileName& subDir,
    const word& cloudName,
    IOstream::streamFormat format
);

//- write generalized field components
template<class Type>
void ensightVolField
(
    const ensightParts& partsList,
    const IOobject& fieldObject,
    const fvMesh& mesh,
    const fileName& dataDir,
    const fileName& subDir,
    IOstream::streamFormat format
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace Foam

#ifdef NoRepository
#   include "ensightOutputFunctions.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
