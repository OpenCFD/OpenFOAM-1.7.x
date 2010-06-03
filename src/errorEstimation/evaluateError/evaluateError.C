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

#include "evaluateError.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "refineCell.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
Foam::evaluateError::evaluateError()
:
    unsplitFaces_(),
    refCells_()
{}


// Construct from components
Foam::evaluateError::evaluateError
(
    const volScalarField& cellError,
    const volVectorField& gradTheta,
    const surfaceScalarField& faceError,
    const labelList& candidateFaces
)
:
    unsplitFaces_(candidateFaces.size()),
    refCells_()
{
    const polyMesh& mesh = cellError.mesh();

    // picks up the error field and the gradient of the variable
    // and appends lists of cells to refine/unrefine based on the width of
    // standard deviation of the error distribution

    // calculate the average error
    scalar avgError = cellError.average().value();

    scalar squareError = sqr(cellError)().average().value();
    scalar deviation = sqrt(squareError - sqr(avgError));

    Info<< "avgError:" << avgError
        << "  squareError:" << squareError
        << "  deviation:" << deviation
        << endl;

    scalar ref = avgError + deviation;
    scalar unref = avgError - deviation;

    Info<< "evaluateError : refinement criterion : " << ref << endl
        << "                unrefinement criterion : " << unref << endl;

    // Coarsen mesh first.
    // Find out set of candidateFaces where error is above crit.

    // Construct to filter unrefinement pattern
//    removeFaces faceRemover(mesh);

    // Keep track of unrefinement pattern.
    boolList markedFace(mesh.nFaces(), false);

    label unsplitFaceI = 0;

    // Subset candidate faces and update refinement pattern interference pattern
    forAll(candidateFaces, candidateFaceI)
    {
        label faceI = candidateFaces[candidateFaceI];

        if (markedFace[faceI])
        {
            Info<< "evaluateError : protected candidate face:" << faceI
                << endl;
        }
        else
        {
//            if (faceError[faceI] < unref)
            if (unsplitFaceI < (candidateFaces.size()/2 + 1))
            {
                unsplitFaces_[unsplitFaceI++] = faceI;

//                faceRemover.markAffectedFaces(faceI, markedFace);
            }
        }
    }

    unsplitFaces_.setSize(unsplitFaceI);

    // Now we have:
    // -unsplitFaces_: all the faces that will be removed 
    // -markedFace   : all the faces affected by this removal.
    // From markedFace protect the cells using them.

    boolList markedCells(mesh.nCells(), false);

//    forAll(markedFace, faceI)
//    {
//        if (markedFace[faceI])
//        {
//            markedCells[mesh.faceOwner()[faceI]] = true;
//
//            if (mesh.isInternalFace(faceI))
//            {
//                markedCells[mesh.faceNeighbour()[faceI]] = true;
//            }
//        }
//    }

    // Select the cells that need to be split.
    // Two pass: count first, select later.

    label refCellI = 0;

    forAll(cellError, cellI)
    {
        if ((cellError[cellI] > ref) && !markedCells[cellI])
        {
            refCellI++;
        }
    }

    refCells_.setSize(refCellI);

    refCellI = 0;

    forAll(cellError, cellI)
    {
        if ((cellError[cellI] > ref) && !markedCells[cellI])
        {
            refCells_[refCellI++] = refineCell(cellI, gradTheta[cellI]);
        }
    }

    Info<< "evaluateError : selected " << unsplitFaces_.size()
        << " faces out of " << candidateFaces.size() << " for removal" << endl;
    Info<< "evaluateError : selected " << refCells_.size()
        << " cells out of " << cellError.size() << " for refinement" << endl;
}


// ************************************************************************* //
