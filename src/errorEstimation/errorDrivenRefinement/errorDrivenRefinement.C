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

\*---------------------------------------------------------------------------*/

#include "errorDrivenRefinement.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "evaluateError.H"
#include "fvc.H"
#include "mapPolyMesh.H"
#include "topoCellLooper.H"
#include "cellCuts.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(errorDrivenRefinement, 0);
    addToRunTimeSelectionTable
    (
        polyMeshModifier,
        errorDrivenRefinement,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::errorDrivenRefinement::errorDrivenRefinement
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyTopoChanger& mme
)
:
    polyMeshModifier(name, index, mme, false),
    refinementEngine_(topoChanger().mesh(), true),
    errorField_(dict.lookup("errorField"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::errorDrivenRefinement::~errorDrivenRefinement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::errorDrivenRefinement::changeTopology() const
{
    const Time& runTime = topoChanger().mesh().time();

    if (runTime.foundObject<volVectorField>(errorField_))
    {
        if (debug)
        {
            Info<< "errorDrivenRefinement::changeTopology() : triggering topo"
                << " change since found errorField "
                << errorField_ << endl;
        }

        return true;
    }
    else
    {
        if (debug)
        {
            Info<< "errorDrivenRefinement::changeTopology() : no topo"
                << " change request from me since no errorField "
                << errorField_ << endl;
        }

        return false;
    }
}


void Foam::errorDrivenRefinement::setRefinement(polyTopoChange& ref) const
{
    // Insert the coarsen/refinement instructions into the topological change

    if (debug)
    {
        Info<< "errorDrivenRefinement::setRefinement(polyTopoChange& ref)"
            << endl;
    }

    const polyMesh& mesh = topoChanger().mesh();

    const Time& runTime = mesh.time();

    if (debug)
    {
        Info<< "Looking up vector field with name " << errorField_ << endl;
    }
    const volVectorField& resError =
        runTime.lookupObject<volVectorField>(errorField_);

    const volScalarField magResError = Foam::mag(resError);

    scalar min = Foam::min(magResError).value();
    scalar max = Foam::max(magResError).value();
    scalar avg = Foam::average(magResError).value();

    if (debug) 
    {
        Info<< "Writing magResError" << endl;
        magResError.write();

        Info<< "min:" << min << " max:" << max << " avg:" << avg << endl;
    }

    // Get faces to remove and cells to refine based on error
    evaluateError refPattern
    (
        magResError,                        // Error on cells
        resError,                           // Error vector on cells
        fvc::interpolate(magResError),      // Error on faces
        refinementEngine_.getSplitFaces()   // Current live split faces
    );


    // Insert mesh refinement into polyTopoChange:
    // - remove split faces
    // - refine cells

    // Give 'hint' of faces to remove to cell splitter.
    const labelList& candidates = refPattern.unsplitFaces();
    ////Hack:no unsplitting
    //labelList candidates;

    labelList removedFaces(refinementEngine_.removeSplitFaces(candidates, ref));

    // Now success will be for every candidates whether face has been removed.
    // Protect cells using face from refinement.

    // List of protected cells
    boolList markedCell(mesh.nCells(), false);

    forAll(removedFaces, i)
    {
        label faceI = removedFaces[i];

        markedCell[mesh.faceOwner()[faceI]] = true;

        if (mesh.isInternalFace(faceI))
        {
            markedCell[mesh.faceNeighbour()[faceI]] = true;
        }
    }
    
    // Repack list of cells to refine.
    List<refineCell> refCells = refPattern.refCells();

    label newRefCellI = 0;

    forAll(refCells, refCellI)
    {
        label cellI = refCells[refCellI].cellNo();

        if (!markedCell[cellI] && (newRefCellI != refCellI))
        {
            refCells[newRefCellI++] = refCells[refCellI];
        }
    }

    if (debug)
    {
        Info<< "errorDrivenRefinement : shrinking refCells from "
            << refCells.size()
            << " to " << newRefCellI << endl;
    }

    refCells.setSize(newRefCellI);

    // Determine cut pattern using topological cell walker
    topoCellLooper cellWalker(mesh);

    cellCuts cuts(mesh, cellWalker, refCells);

    // Do actual splitting
    refinementEngine_.setRefinement(cuts, ref);
}


// Has the responsability of moving my newly introduced points onto the right
// place. This is since the whole mesh might e.g. have been moved by another
// meshmodifier. So using preMotionPoints is hack for if I am only meshModifier.
// Good solution:
// - remember new point label of introduced point and vertices
// of edge it is created from (in setRefinement)
// - in here reposition point at correct position between current vertex
// position of edge endpoints.
void Foam::errorDrivenRefinement::modifyMotionPoints
(
    pointField& motionPoints
) const
{
    if (debug)
    {
        Info<< "errorDrivenRefinement::modifyMotionPoints(*pointField&)" << endl;
    }
}


void Foam::errorDrivenRefinement::updateMesh(const mapPolyMesh& morphMap)
{
    // Mesh has changed topologically. Update local topological data
    if (debug)
    {
        Info<< "errorDrivenRefinement::updateMesh"
            << "(const mapPolyMesh& morphMap)" << endl;
    }
    refinementEngine_.updateMesh(morphMap);
}


void Foam::errorDrivenRefinement::write(Ostream& os) const
{
    os  << nl << type() << nl;
}


void Foam::errorDrivenRefinement::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl
        << "    type " << type()
        << token::END_STATEMENT << nl
        << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
