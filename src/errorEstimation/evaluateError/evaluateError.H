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

Class
    Foam::evaluateError

Description
    Foam::evaluateError

SourceFiles
    evaluateError.C

\*---------------------------------------------------------------------------*/

#ifndef evaluateError_H
#define evaluateError_H

#include "labelList.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class refineCell;

/*---------------------------------------------------------------------------*\
                           Class evaluateError Declaration
\*---------------------------------------------------------------------------*/

class evaluateError
{
    // Private data

        // splitFaces to remove.
        labelList unsplitFaces_;

        // cells to refine.
        List<refineCell> refCells_;


    // Private Member Functions


public:

    // Constructors

        //- Construct null
        evaluateError();


        //- Construct from error value, error direction. Determine faces to
        //  remove (out of candidate faces) and cells to refine (all cells):
        //
        //  remove all faces where the (face)error is less than
        //     average - standard deviation
        //  refine all cells where the (cell)error is more than
        //      average + standard deviation. gradTheta is used to determine
        //      direction to refine in.
        evaluateError
        (
            const volScalarField& cellError,
            const volVectorField& gradTheta,
            const surfaceScalarField& faceError,
            const labelList& candidateFaces
        );


    // Member Functions

        const labelList& unsplitFaces() const
        {
            return unsplitFaces_;
        }

        const List<refineCell>& refCells() const
        {
            return refCells_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
