/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    Foam::patchInteractionData

Description
    Helper class for the LocalInteraction patch interaction model

\*---------------------------------------------------------------------------*/

#ifndef patchInteractionData_H
#define patchInteractionData_H

#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class patchInteractionData Declaration
\*---------------------------------------------------------------------------*/

// Forward declaration of classes
class patchInteractionData;

// Forward declaration of friend functions
Istream& operator>>
(
    Istream& is,
    patchInteractionData& pid
);


class patchInteractionData
{
    // Private data

        //- Interaction type name
        word interactionTypeName_;

        //- Patch name
        word patchName_;

        //- Elasticity coefficient
        scalar e_;

        //- Restitution coefficient
        scalar mu_;


public:

    // Constructor

        //- Construct null
        patchInteractionData();


    // Member functions

        // Access

            //- Return const access to the interaction type name
            const word& interactionTypeName() const;

            //- Return const access to the patch name
            const word& patchName() const;

            //- Return const access to the elasticity coefficient
            scalar e() const;

            //- Return const access to the restitution coefficient
            scalar mu() const;


        // I-O

            //- Istream operator
            friend Istream& operator>>
            (
                Istream& is,
                patchInteractionData& pid
            );
/*            {
                is.check
                (
                    "Istream& operator>>"
                    "(Istream&, patchInteractionData&)"
                );

                const dictionaryEntry entry(dictionary::null, is);

                pid.patchName_ = entry.keyword();
                entry.lookup("type") >> pid.interactionTypeName_;
                pid.e_ = entry.lookupOrDefault<scalar>("e", 1.0);
                pid.mu_ = entry.lookupOrDefault<scalar>("mu", 0.0);

                if
                (
                    PatchInteractionModel<CloudType>::wordToInteractionType
                    (
                        pid.interactionTypeName_
                    )
                == PatchInteractionModel<CloudType>::itOther)
                {
                    FatalErrorIn
                    (
                        "friend Istream& operator>>"
                        "("
                            "Istream&, "
                            "patchInteractionData&"
                        ")"
                    )   << "Unknown patch interaction type "
                        << pid.interactionTypeName_
                        << ". Valid selections are:"
                        << PatchInteractionModel<CloudType>::
                        interactionTypeNames_
                        << endl << abort(FatalError);
                }

                return is;
            }
*/};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
