/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::referredMolecule

Description

SourceFiles
    referredMoleculeI.H
    referredMolecule.C

\*---------------------------------------------------------------------------*/

#ifndef referredMolecule_H
#define referredMolecule_H

#include "vector.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class referredMolecule Declaration
\*---------------------------------------------------------------------------*/

class referredMolecule
{
    // Private data

        label id_;

        vector position_;

        List<vector> sitePositions_;


public:

    // Constructors

        //- Construct null
        referredMolecule();

        //- Construct from components
        referredMolecule
        (
            const label id,
            const vector& position,
            const List<vector>& sitePositions
        );


    // Destructor

        virtual ~referredMolecule();


    // Member Functions

        // Access

            inline label id() const;

            inline const vector& position() const;

            inline const List<vector>& sitePositions() const;


    // Friend Operators

        inline friend bool operator==
        (
            const referredMolecule& a,
            const referredMolecule& b
        );

        inline friend bool operator!=
        (
            const referredMolecule& a,
            const referredMolecule& b
        );


    // IOstream Operators

        friend Istream& operator>>
        (
            Istream&,
            referredMolecule&
        );

        friend Ostream& operator<<
        (
            Ostream&,
            const referredMolecule&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "referredMoleculeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
