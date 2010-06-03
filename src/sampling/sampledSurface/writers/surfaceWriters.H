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

InClass
    Foam::surfaceWriters

Description

\*---------------------------------------------------------------------------*/

#ifndef surfaceWriters_H
#define surfaceWriters_H

#include "surfaceWriter.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Only used internally
#define makeTypeSurfaceWritersTypeName(typeWriter, dataType)                 \
                                                                             \
    defineNamedTemplateTypeNameAndDebug(typeWriter< dataType >, 0)

// Sometimes used externally
#define makeSurfaceWritersTypeName(typeWriter)                               \
                                                                             \
    makeTypeSurfaceWritersTypeName(typeWriter, scalar);                      \
    makeTypeSurfaceWritersTypeName(typeWriter, vector);                      \
    makeTypeSurfaceWritersTypeName(typeWriter, sphericalTensor);             \
    makeTypeSurfaceWritersTypeName(typeWriter, symmTensor);                  \
    makeTypeSurfaceWritersTypeName(typeWriter, tensor)

// Define type info for single dataType template instantiation (eg, vector)
#define makeSurfaceWriterType(typeWriter, dataType)                          \
                                                                             \
    defineNamedTemplateTypeNameAndDebug(typeWriter< dataType >, 0);          \
    addTemplatedToRunTimeSelectionTable                                      \
    (                                                                        \
        surfaceWriter, typeWriter, dataType, word                            \
    )


// Define type info for scalar, vector etc. instantiations
#define makeSurfaceWriters(typeWriter)                                       \
                                                                             \
    makeSurfaceWriterType(typeWriter, scalar);                               \
    makeSurfaceWriterType(typeWriter, vector);                               \
    makeSurfaceWriterType(typeWriter, sphericalTensor);                      \
    makeSurfaceWriterType(typeWriter, symmTensor);                           \
    makeSurfaceWriterType(typeWriter, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
