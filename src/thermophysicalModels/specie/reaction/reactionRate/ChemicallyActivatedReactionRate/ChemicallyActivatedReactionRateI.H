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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class ReactionRate, class ChemicallyActivationFunction>
inline
ChemicallyActivatedReactionRate<ReactionRate, ChemicallyActivationFunction>::
ChemicallyActivatedReactionRate
(
    const ReactionRate& k0,
    const ReactionRate& kInf,
    const ChemicallyActivationFunction& F,
    const thirdBodyEfficiencies& tbes
)
:
    k0_(k0),
    kInf_(kInf),
    F_(F),
    thirdBodyEfficiencies_(tbes)
{}


//- Construct from Istream
template<class ReactionRate, class ChemicallyActivationFunction>
inline
ChemicallyActivatedReactionRate<ReactionRate, ChemicallyActivationFunction>::
ChemicallyActivatedReactionRate
(
    const speciesTable& species,
    Istream& is
)
:
    k0_(species, is.readBegin("ChemicallyActivatedReactionRate(Istream&)")),
    kInf_(species, is),
    F_(is),
    thirdBodyEfficiencies_(species, is)
{
    is.readEnd("ChemicallyActivatedReactionRate(Istream&)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionRate, class ChemicallyActivationFunction>
inline scalar ChemicallyActivatedReactionRate
    <ReactionRate, ChemicallyActivationFunction>::operator()
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    scalar k0 = k0_(T, p, c);
    scalar kInf = kInf_(T, p, c);

    scalar Pr = k0*thirdBodyEfficiencies_.M(c)/kInf;

    return k0*(1/(1 + Pr))*F_(T, Pr);
}


template<class ReactionRate, class ChemicallyActivationFunction>
inline Ostream& operator<<
(
    Ostream& os,
    const ChemicallyActivatedReactionRate
        <ReactionRate, ChemicallyActivationFunction>& carr
)
{
    os  << token::BEGIN_LIST
        << carr.k0_ << token::SPACE << carr.kInf_ << token::SPACE << carr.F_
        << token::END_LIST;
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
