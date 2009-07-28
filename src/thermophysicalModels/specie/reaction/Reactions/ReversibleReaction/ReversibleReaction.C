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

Description
    
\*---------------------------------------------------------------------------*/

#include "ReversibleReaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class ReactionThermo, class ReactionRate>
ReversibleReaction<ReactionThermo, ReactionRate>::ReversibleReaction
(
    const Reaction<ReactionThermo>& reaction,
    const ReactionRate& k
)
:
    Reaction<ReactionThermo>(reaction),
    k_(k)
{}


// Construct from components
template<class ReactionThermo, class ReactionRate>
ReversibleReaction<ReactionThermo, ReactionRate>::ReversibleReaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    Istream& is
)
:
    Reaction<ReactionThermo>(species, thermoDatabase, is),
    k_(species, is)
{}


// Construct as copy given new speciesTable
template<class ReactionThermo, class ReactionRate>
ReversibleReaction<ReactionThermo, ReactionRate>::ReversibleReaction
(
    const ReversibleReaction<ReactionThermo, ReactionRate>& rr,
    const speciesTable& species
)
:
    Reaction<ReactionThermo>(rr, species),
    k_(rr.k_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo, class ReactionRate>
scalar ReversibleReaction<ReactionThermo, ReactionRate>::kf
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return k_(T, p, c);
}


template<class ReactionThermo, class ReactionRate>
scalar ReversibleReaction<ReactionThermo, ReactionRate>::kr
(
    const scalar kfwd,
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return kfwd/this->Kc(T);
}


template<class ReactionThermo, class ReactionRate>
scalar ReversibleReaction<ReactionThermo, ReactionRate>::kr
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return kr(kf(T, p, c), T, p, c);
}


template<class ReactionThermo, class ReactionRate>
void ReversibleReaction<ReactionThermo, ReactionRate>::write
(
    Ostream& os
) const
{
    Reaction<ReactionThermo>::write(os);
    os  << token::SPACE << k_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
