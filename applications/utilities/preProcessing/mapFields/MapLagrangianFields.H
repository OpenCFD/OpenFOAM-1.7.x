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

InNamespace
    Foam

Description
    Gets the indices of (source)particles that have been appended to the
    target cloud and maps the lagrangian fields accordingly.

\*---------------------------------------------------------------------------*/

#ifndef MapLagrangianFields_H
#define MapLagrangianFields_H

#include "cloud.H"
#include "GeometricField.H"
#include "meshToMesh.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Gets the indices of (source)particles that have been appended to the
//  target cloud and maps the lagrangian fields accordingly.
template<class Type>
void MapLagrangianFields
(
    const string& cloudName,
    const IOobjectList& objects,
    const meshToMesh& meshToMeshInterp,
    const labelList& addParticles
)
{
    const fvMesh& meshTarget = meshToMeshInterp.toMesh();

    IOobjectList fields = objects.lookupClass(IOField<Type>::typeName);

    forAllIter(IOobjectList, fields, fieldIter)
    {
        Info<< "    mapping lagrangian field " << fieldIter()->name() << endl;

        // Read field (does not need mesh)
        IOField<Type> fieldSource(*fieldIter());

        // Map
        IOField<Type> fieldTarget
        (
            IOobject
            (
                fieldIter()->name(),
                meshTarget.time().timeName(),
                cloud::prefix/cloudName,
                meshTarget,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            addParticles.size()
        );
        forAll(addParticles, i)
        {
            fieldTarget[i] = fieldSource[addParticles[i]];
        }

        // Write field
        fieldTarget.write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
