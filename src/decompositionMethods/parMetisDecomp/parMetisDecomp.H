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

Class
    Foam::parMetisDecomp

Description

SourceFiles
    parMetisDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef parMetisDecomp_H
#define parMetisDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class parMetisDecomp Declaration
\*---------------------------------------------------------------------------*/

class parMetisDecomp
:
    public decompositionMethod
{
    // Private data

        const polyMesh& mesh_;


    // Private Member Functions

        //- Insert list in front of list.
        template<class Type>
        static void prepend(const UList<Type>&, List<Type>&);
        //- Insert list at end of list.
        template<class Type>
        static void append(const UList<Type>&, List<Type>&);

        label decompose
        (
            Field<int>& xadj,
            Field<int>& adjncy,
            const pointField& cellCentres,
            Field<int>& cellWeights,
            Field<int>& faceWeights,
            const List<int>& options,

            List<int>& finalDecomp
        );


        //- Disallow default bitwise copy construct and assignment
        void operator=(const parMetisDecomp&);
        parMetisDecomp(const parMetisDecomp&);


public:

    //- Runtime type information
    TypeName("parMetis");


    // Constructors

        //- Construct given the decomposition dictionary and mesh
        parMetisDecomp
        (
            const dictionary& decompositionDict,
            const polyMesh& mesh
        );


    // Destructor

        ~parMetisDecomp()
        {}


    // Member Functions

        //- parMetis handles Foam processor boundaries
        virtual bool parallelAware() const
        {
            return true;
        }

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        virtual labelList decompose(const pointField&);

        //- Return for every coordinate the wanted processor number. Gets
        //  passed agglomeration map (from fine to coarse cells) and coarse cell
        //  location. Can be overridden by decomposers that provide this
        //  functionality natively.
        virtual labelList decompose
        (
            const labelList& agglom,
            const pointField&
        );

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided mesh connectivity.
        //  The connectivity is equal to mesh.cellCells() except for
        //  - in parallel the cell numbers are global cell numbers (starting
        //    from 0 at processor0 and then incrementing all through the
        //    processors)
        //  - the connections are across coupled patches
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "parMetisDecompTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
