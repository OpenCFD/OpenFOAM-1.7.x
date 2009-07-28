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
    Foam::polyPatchID

Description
    A class holds the data needed to identify a patch in a dynamic mesh.

    The patch is identified by name and its index in the boundary mesh
    is updated if the mesh has changed.

\*---------------------------------------------------------------------------*/

#ifndef polyPatchID_H
#define polyPatchID_H

#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class polyPatchID;
Ostream& operator<<(Ostream& os, const polyPatchID& p);


/*---------------------------------------------------------------------------*\
                           Class polyPatchID Declaration
\*---------------------------------------------------------------------------*/

class polyPatchID
{
    // Private data

        //- Patch name
        word name_;

        //- Patch index
        label index_;


public:

    // Constructors

        //- Construct from name
        polyPatchID(const word& name, const polyBoundaryMesh& bm)
        :
            name_(name),
            index_(bm.findPatchID(name))
        {}

        //- Construct from Istream
        polyPatchID(Istream& is, const polyBoundaryMesh& bm)
        :
            name_(is),
            index_(bm.findPatchID(name_))
        {}


    // Member Functions

        // Access

            //- Return name
            const word& name() const
            {
                return name_;
            }

            //- Return index
            label index() const
            {
                return index_;
            }

            //- Has the patch been found
            bool active() const
            {
                return index_ > -1;
            }


        // Edit

            //- Update
            void update(const polyBoundaryMesh& bm)
            {
                index_ = bm.findPatchID(name_);
            }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const polyPatchID& p)
        {
            os  << token::BEGIN_LIST
                << p.name_ << token::SPACE
                << p.index_
                << token::END_LIST;

            // Check state of Ostream
            os.check("Ostream& operator<<(Ostream&, const polyPatchID&)");

            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
