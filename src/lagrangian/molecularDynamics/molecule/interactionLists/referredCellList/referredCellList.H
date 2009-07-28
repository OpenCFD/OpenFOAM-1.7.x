/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    Foam::referredCellList

Description

SourceFiles
    referredCellListI.H
    referredCellList.C

\*---------------------------------------------------------------------------*/

#ifndef referredCellList_H
#define referredCellList_H

#include "referredCell.H"
#include "molecule.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class interactionLists;

/*---------------------------------------------------------------------------*\
                      Class referredCellList Declaration
\*---------------------------------------------------------------------------*/

class referredCellList
:
    public List<referredCell>
{
    // Private data

        const interactionLists& il_;


    // Private Member Functions

        void buildReferredCellList
        (
            bool pointPointListBuild
        );


public:

    // Constructors

        //- Construct lists by searching the mesh
        referredCellList
        (
            interactionLists& il,
            bool pointPointListBuild
        );

        //- Construct from file
        referredCellList (interactionLists& il);


    // Destructor

        ~referredCellList();


    // Member Functions

        void referMolecules(const List<DynamicList<molecule*> >& cellOccupancy);

        inline const interactionLists& il() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "referredCellListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
