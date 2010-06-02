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
    Foam::simpleGeomDecomp

Description

SourceFiles
    simpleGeomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef simpleGeomDecomp_H
#define simpleGeomDecomp_H

#include "geomDecomp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class simpleGeomDecomp Declaration
\*---------------------------------------------------------------------------*/

class simpleGeomDecomp
:
    public geomDecomp
{
    // Private Member Functions

        void assignToProcessorGroup(labelList& processorGroup, const label);

        void assignToProcessorGroup
        (
            labelList& processorGroup,
            const label nProcGroup,
            const labelList& indices,
            const scalarField& weights,
            const scalar summedWeights
        );

        //- Disallow default bitwise copy construct and assignment
        void operator=(const simpleGeomDecomp&);
        simpleGeomDecomp(const simpleGeomDecomp&);


public:

    //- Runtime type information
    TypeName("simple");


    // Constructors

        //- Construct given the decomposition dictionary
        simpleGeomDecomp(const dictionary& decompositionDict);

        //- Construct given the decomposition dictionary and mesh
        simpleGeomDecomp
        (
            const dictionary& decompositionDict,
            const polyMesh& mesh
        );


    // Destructor

        virtual ~simpleGeomDecomp()
        {}


    // Member Functions

        virtual bool parallelAware() const
        {
            // simpleDecomp does not attempt to do anything across proc
            // boundaries
            return false;
        }

        virtual labelList decompose
        (
            const pointField& points
        );

        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& pointWeights
        );

        //- Explicitly provided connectivity
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        )
        {
            return decompose(cc, cWeights);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
