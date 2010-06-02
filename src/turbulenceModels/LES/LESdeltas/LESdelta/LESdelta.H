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
    Foam::LESdelta

Description
    Abstract base class for LES deltas

SourceFiles
    LESdelta.C
    newDelta.C

\*---------------------------------------------------------------------------*/

#ifndef LESdelta_H
#define LESdelta_H

#include "volFields.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                          Class LESdelta Declaration
\*---------------------------------------------------------------------------*/

class LESdelta
{

protected:

    // Protected data

        const fvMesh& mesh_;

        volScalarField delta_;


private:

    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        LESdelta(const LESdelta&);
        void operator=(const LESdelta&);


public:

    //- Runtime type information
    TypeName("LESdelta");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            LESdelta,
            dictionary,
            (
                const word& name,
                const fvMesh& mesh,
                const dictionary& LESdeltaDict
            ),
            (name, mesh, LESdeltaDict)
        );


    // Constructors

        //- Construct from name and mesh
        LESdelta(const word& name, const fvMesh&);


    // Selectors

        //- Return a reference to the selected LES delta
        static autoPtr<LESdelta> New
        (
            const word& name,
            const fvMesh&,
            const dictionary&
        );


    //- Destructor
    virtual ~LESdelta()
    {}


    // Member Functions

        //- Return mesh reference
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Read the LESdelta dictionary
        virtual void read(const dictionary&) = 0;

        // Correct values
        virtual void correct() = 0;


    // Member Operators

        virtual operator const volScalarField&() const
        {
            return delta_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
