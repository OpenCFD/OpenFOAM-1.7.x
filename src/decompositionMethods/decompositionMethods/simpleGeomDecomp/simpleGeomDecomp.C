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


void Foam::simpleGeomDecomp::assignToProcessorGroup
(
    labelList& processorGroup,
    const label nProcGroup,
    const labelList& indices,
    const scalarField& weights,
    const scalar summedWeights
)
{
    // This routine gets the sorted points.
    // Easiest to explain with an example.
    // E.g. 400 points, sum of weights : 513.
    // Now with number of divisions in this direction (nProcGroup) : 4
    // gives the split at 513/4 = 128
    // So summed weight from 0..128 goes into bin 0,
    //     ,,              128..256 goes into bin 1
    //   etc.
    // Finally any remaining ones go into the last bin (3).

    const scalar jump = summedWeights/nProcGroup;
    const label nProcGroupM1 = nProcGroup - 1;
    scalar sumWeights = 0;
    label ind = 0;
    label j = 0;

    // assign cells to all except last group.
    for (j=0; j<nProcGroupM1; j++)
    {
        const scalar limit = jump*scalar(j + 1);
        while (sumWeights < limit)
        {
            sumWeights += weights[indices[ind]];
            processorGroup[ind++] = j;
        }
    }
    // Ensure last included.
    while (ind < processorGroup.size())
    {
       processorGroup[ind++] = nProcGroupM1;
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
    forAll(pointIndices, i)
    {
        pointIndices[i] = i;
    }

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


Foam::labelList Foam::simpleGeomDecomp::decompose
(
    const pointField& points,
    const scalarField& weights
)
{
    // construct a list for the final result
    labelList finalDecomp(points.size());

    labelList processorGroups(points.size());

    labelList pointIndices(points.size());
    forAll(pointIndices, i)
    {
        pointIndices[i] = i;
    }

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

    const scalar summedWeights = sum(weights);
    assignToProcessorGroup
    (
        processorGroups,
        n_.x(),
        pointIndices,
        weights,
        summedWeights
    );

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

    assignToProcessorGroup
    (
        processorGroups,
        n_.y(),
        pointIndices,
        weights,
        summedWeights
    );

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

    assignToProcessorGroup
    (
        processorGroups,
        n_.z(),
        pointIndices,
        weights,
        summedWeights
    );

    forAll(points, i)
    {
        finalDecomp[pointIndices[i]] += n_.x()*n_.y()*processorGroups[i];
    }

    return finalDecomp;
}
// ************************************************************************* //
