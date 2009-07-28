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

Description

\*---------------------------------------------------------------------------*/

#include "simpleGeomDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleGeomDecomp, 0);

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        simpleGeomDecomp,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        decompositionMethod,
        simpleGeomDecomp,
        dictionaryMesh
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// assignToProcessorGroup : given nCells cells and nProcGroup processor
// groups to share them, how do we share them out? Answer : each group
// gets nCells/nProcGroup cells, and the first few get one
// extra to make up the numbers. This should produce almost
// perfect load balancing

void Foam::simpleGeomDecomp::assignToProcessorGroup
(
    labelList& processorGroup,
    const label nProcGroup
)
{
    label jump = processorGroup.size()/nProcGroup;
    label jumpb = jump + 1;
    label fstProcessorGroup = processorGroup.size() - jump*nProcGroup;

    label ind = 0;
    label j = 0;

    // assign cells to the first few processor groups (those with
    // one extra cell each
    for (j=0; j<fstProcessorGroup; j++)
    {
        for (register label k=0; k<jumpb; k++)
        {
            processorGroup[ind++] = j;
        }
    }

    // and now to the `normal' processor groups
    for (; j<nProcGroup; j++)
    {
        for (register label k=0; k<jump; k++)
        {
            processorGroup[ind++] = j;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleGeomDecomp::simpleGeomDecomp(const dictionary& decompositionDict)
:
    geomDecomp(decompositionDict, typeName)
{}


Foam::simpleGeomDecomp::simpleGeomDecomp
(
    const dictionary& decompositionDict,
    const polyMesh&
)
:
    geomDecomp(decompositionDict, typeName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::simpleGeomDecomp::decompose(const pointField& points)
{
    // construct a list for the final result
    labelList finalDecomp(points.size());

    labelList processorGroups(points.size());

    labelList pointIndices(points.size());
    forAll(pointIndices, i) {pointIndices[i] = i;}

    pointField rotatedPoints = rotDelta_ & points;

    // and one to take the processor group id's. For each direction.
    // we assign the processors to groups of processors labelled
    // 0..nX to give a banded structure on the mesh. Then we
    // construct the actual processor number by treating this as
    // the units part of the processor number.
    sort
    (
        pointIndices,
        UList<scalar>::less(rotatedPoints.component(vector::X))
    );

    assignToProcessorGroup(processorGroups, n_.x());

    forAll(points, i)
    {
        finalDecomp[pointIndices[i]] = processorGroups[i];
    }


    // now do the same thing in the Y direction. These processor group
    // numbers add multiples of nX to the proc. number (columns)
    sort
    (
        pointIndices,
        UList<scalar>::less(rotatedPoints.component(vector::Y))
    );

    assignToProcessorGroup(processorGroups, n_.y());

    forAll(points, i)
    {
        finalDecomp[pointIndices[i]] += n_.x()*processorGroups[i];
    }


    // finally in the Z direction. Now we add multiples of nX*nY to give
    // layers
    sort
    (
        pointIndices,
        UList<scalar>::less(rotatedPoints.component(vector::Z))
    );

    assignToProcessorGroup(processorGroups, n_.z());

    forAll(points, i)
    {
        finalDecomp[pointIndices[i]] += n_.x()*n_.y()*processorGroups[i];
    }


    return finalDecomp;
}


// ************************************************************************* //
