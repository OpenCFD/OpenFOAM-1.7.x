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

#include "foamFileSurfaceWriter.H"

#include "OFstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::foamFileSurfaceWriter<Type>::foamFileSurfaceWriter()
:
    surfaceWriter<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::foamFileSurfaceWriter<Type>::~foamFileSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::foamFileSurfaceWriter<Type>::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const bool verbose
) const
{
    fileName surfaceDir(outputDir/surfaceName);

    if (!isDir(surfaceDir))
    {
        mkDir(surfaceDir);
    }

    if (verbose)
    {
        Info<< "Writing geometry to " << surfaceDir << endl;
    }

    // Points
    OFstream(surfaceDir/"points")() << points;

    // Faces
    OFstream(surfaceDir/"faces")() << faces;
}


template<class Type>
void Foam::foamFileSurfaceWriter<Type>::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const fileName& fieldName,
    const Field<Type>& values,
    const bool verbose
) const
{
    fileName surfaceDir(outputDir/surfaceName);

    if (!isDir(surfaceDir))
    {
        mkDir(surfaceDir);
    }

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << surfaceDir << endl;
    }

    // geometry should already have been written

    // Values to separate directory (e.g. "scalarField/p")

    fileName foamName(pTraits<Type>::typeName);
    fileName valuesDir(surfaceDir  / (foamName + Field<Type>::typeName));

    if (!isDir(valuesDir))
    {
        mkDir(valuesDir);
    }

    // values
    OFstream(valuesDir/fieldName)() << values;
}


// ************************************************************************* //
