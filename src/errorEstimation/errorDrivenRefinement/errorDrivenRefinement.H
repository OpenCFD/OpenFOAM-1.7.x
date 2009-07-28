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
    Foam::errorDrivenRefinement

Description
    Refines and coarsens based on error estimate.

SourceFiles
    errorDrivenRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef errorDrivenRefinement_H
#define errorDrivenRefinement_H

#include "polyMeshModifier.H"
#include "undoableMeshCutter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                       Class errorDrivenRefinement Declaration
\*---------------------------------------------------------------------------*/

class errorDrivenRefinement
:
    public polyMeshModifier
{
    // Private data

        //- Refinement/coarsening engine
        mutable undoableMeshCutter refinementEngine_;

        //- Name of volVectorField which contains error.
        word errorField_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        errorDrivenRefinement(const errorDrivenRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const errorDrivenRefinement&);


public:

    //- Runtime type information
    TypeName("errorDrivenRefinement");


    // Constructors

        //- Construct from dictionary
        errorDrivenRefinement
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyTopoChanger& mme
        );


    // Destructor

        virtual ~errorDrivenRefinement();


    // Member Functions

        //- Check for topology change
        virtual bool changeTopology() const;

        //- Insert the layer addition/removal instructions
        //  into the topological change
        virtual void setRefinement(polyTopoChange&) const;

        //- Modify motion points to comply with the topological change
        virtual void modifyMotionPoints(pointField& motionPoints) const;

        //- Force recalculation of locally stored data on topological change
        virtual void updateMesh(const mapPolyMesh&);

        //- Write
        virtual void write(Ostream&) const;

        //- Write dictionary
        virtual void writeDict(Ostream&) const;


    // Access

        //- Underlying mesh modifier        
        const undoableMeshCutter& refinementEngine() const
        {
            return refinementEngine_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
