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
    Foam::treeElem

Description
    Common functionality of treeNode and treeLeaf.

SourceFiles
    treeElem.C

\*---------------------------------------------------------------------------*/

#ifndef treeElem_H
#define treeElem_H

#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class treeElem Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class treeElem
{
    // Private data

        //- Bounding box of this node
        treeBoundBox bb_;

public:

    // Constructors

        //- Construct from bounding box
        treeElem(const treeBoundBox& bb)
        :
            bb_(bb)
        {}


    // Member Functions

        // Access

            //- Bounding box of this node
            const treeBoundBox& bb() const
            {
                return bb_;
            }

            //- Bounding box of this node
            treeBoundBox& bb()
            {
                return bb_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
