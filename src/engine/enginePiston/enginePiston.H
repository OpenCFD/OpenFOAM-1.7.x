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
    Foam::enginePiston

Description
    Foam::enginePiston

SourceFiles
    enginePiston.C

\*---------------------------------------------------------------------------*/

#ifndef enginePiston_H
#define enginePiston_H

#include "polyPatchID.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class engineTime;

/*---------------------------------------------------------------------------*\
                           Class enginePiston Declaration
\*---------------------------------------------------------------------------*/

class enginePiston
{
    // Private data

        //- Reference to engine mesh
        const polyMesh& mesh_;

        //- Reference to engine database
        const engineTime& engineDB_;

        //- Piston patch
        polyPatchID patchID_;

        //- Coordinate system
        autoPtr<coordinateSystem> csPtr_;


        // Piston layering data

            //- Min layer thickness
            const scalar minLayer_;

            //- Max layer thickness
            const scalar maxLayer_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        enginePiston(const enginePiston&);

        //- Disallow default bitwise assignment
        void operator=(const enginePiston&);


public:

    // Static data members


    // Constructors

        //- Construct from components
        enginePiston
        (
            const polyMesh& mesh,
            const word& pistonPatchName,
            const autoPtr<coordinateSystem>& pistonCS,
            const scalar minLayer,
            const scalar maxLayer

        );

        //- Construct from dictionary
        enginePiston
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor - default


    // Member Functions

        //- Return coordinate system
        const coordinateSystem& cs() const
        {
            return csPtr_();
        }

        //- Return ID of piston patch
        const polyPatchID& patchID() const
        {
            return patchID_;
        }

        // Piston layering thickness

            scalar minLayer() const
            {
                return minLayer_;
            }

            scalar maxLayer() const
            {
                return maxLayer_;
            }


        //- Write dictionary
        void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
