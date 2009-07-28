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
    Foam::LocalInteraction

Description
    Patch interaction specified on a patch-by-patch basis

\*---------------------------------------------------------------------------*/

#ifndef LocalInteraction_H
#define LocalInteraction_H

#include "PatchInteractionModel.H"
#include "dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class LocalInteraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class LocalInteraction
:
    public PatchInteractionModel<CloudType>
{
    class patchInteractionData
    {
        // Private data

            //- Patch name
            word patchName_;

            //- Elasticity coefficient
            scalar e_;

            //- Restitution coefficient
            scalar mu_;


    public:

        //- Construct null
        patchInteractionData()
        :
            patchName_("unknownPatch"),
            e_(0.0),
            mu_(0.0)
        {}

        //- Construct from dictionary
        patchInteractionData(const dictionary& dict);

        // Member functions

            // Access

                //- Return const access to the patch name
                const word& patchName() const
                {
                    return patchName_;
                }

                //- Return const access to the elasticity coefficient
                scalar e() const
                {
                    return e_;
                }

                //- Return const access to the restitution coefficient
                scalar mu() const
                {
                    return mu_;
                }


            // I-O

                //- Istream operator
                friend Istream& operator>>(Istream& is, patchInteractionData& pid)
                {
                    is.check
                    (
                        "Istream& operator>>"
                        "(Istream&, patchInteractionData&)"
                    );

                    const dictionaryEntry entry(dictionary::null, is);

                    pid.patchName_ = entry.keyword();
                    entry.lookup("e") >> pid.e_;
                    entry.lookup("mu") >> pid.mu_;

                    return is;
                }
    };


    // Private data

        //- List of participating patches
        const List<patchInteractionData> patchData_;

        //- List of participating patch ids
        List<label> patchIds_;


    // Private member functions

        //- Returns local patchI if patch is in patchIds_ list
        label applyToPatch(const label globalPatchI) const;



public:

    //- Runtime type information
    TypeName("LocalInteraction");


    // Constructors

        //- Construct from dictionary
        LocalInteraction(const dictionary& dict, CloudType& cloud);


    //- Destructor
    virtual ~LocalInteraction();


    // Member Functions

        //- Flag to indicate whether model activates patch interaction model
        virtual bool active() const;

        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            const polyPatch& pp,
            const label faceId,
            vector& U
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "LocalInteraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
