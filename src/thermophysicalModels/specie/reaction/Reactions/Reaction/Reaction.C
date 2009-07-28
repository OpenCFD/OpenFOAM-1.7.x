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
    Simple extension of ReactionThermo to handle Reaction kinetics in addition
    to the equilibrium thermodynamics already handled.

\*---------------------------------------------------------------------------*/

#include "Reaction.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ReactionThermo>
void Reaction<ReactionThermo>::setThermo
(
    const HashPtrTable<ReactionThermo>& thermoDatabase
)
{
    ReactionThermo::operator=
    (
        rhs_[0].stoichCoeff*(*thermoDatabase[species_[rhs_[0].index]])
    );

    for (label i=1; i<rhs_.size(); i++)
    {
        operator+=
        (
            rhs_[i].stoichCoeff*(*thermoDatabase[species_[rhs_[i].index]])
        );
    }

    for (label i=0; i<lhs_.size(); i++)
    {
        operator-=
        (
            lhs_[i].stoichCoeff*(*thermoDatabase[species_[lhs_[i].index]])
        );
    }
}


// Construct from components
template<class ReactionThermo>
Reaction<ReactionThermo>::Reaction
(
    const speciesTable& species,
    const List<specieCoeffs>& lhs,
    const List<specieCoeffs>& rhs,
    const HashPtrTable<ReactionThermo>& thermoDatabase
)
:
    ReactionThermo(*thermoDatabase[species[0]]),
    species_(species),
    lhs_(lhs),
    rhs_(rhs)
{
    setThermo(thermoDatabase);
}


// Construct as copy given new speciesTable
template<class ReactionThermo>
Reaction<ReactionThermo>::Reaction
(
    const Reaction<ReactionThermo>& r,
    const speciesTable& species
)
:
    ReactionThermo(r),
    species_(species),
    lhs_(r.lhs_),
    rhs_(r.rhs_)
{}


template<class ReactionThermo>
Reaction<ReactionThermo>::specieCoeffs::specieCoeffs
(
    const speciesTable& species,
    Istream& is
)
{
    token t(is);

    if (t.isNumber())
    {
        stoichCoeff = t.number();
        is >> t;
    }
    else
    {
        stoichCoeff = 1.0;
    }

    exponent = stoichCoeff;

    if (t.isWord())
    {
        word specieName = t.wordToken();

        size_t i = specieName.find('^');

        if (i != word::npos)
        {
            string exponentStr = specieName
            (
                i + 1,
                specieName.size() - i - 1
            );
            exponent = atof(exponentStr.c_str());
            specieName = specieName(0, i);
        }

        index = species[specieName];
    }
    else
    {
        FatalIOErrorIn("Reaction<ReactionThermo>::lrhs(Istream& is)", is)
            << "Expected a word but found " << t.info()
            << exit(FatalIOError);
    }
}


template<class ReactionThermo>
void Reaction<ReactionThermo>::setLRhs(Istream& is)
{
    DynamicList<specieCoeffs> dlrhs;

    while (is)
    {
        dlrhs.append(specieCoeffs(species_, is));

        token t(is);

        if (t.isPunctuation())
        {
            if (t == token::ADD)
            {
            }
            else if (t == token::ASSIGN)
            {
                lhs_ = dlrhs.shrink();
                dlrhs.clear();
            }
            else
            {
                rhs_ = dlrhs.shrink();
                is.putBack(t);
                return;
            }
        }
        else
        {
            rhs_ = dlrhs.shrink();
            is.putBack(t);
            return;
        }
    }

    FatalIOErrorIn("Reaction<ReactionThermo>::lrhs(Istream& is)", is)
        << "Cannot continue reading reaction data from stream"
        << exit(FatalIOError);
}


//- Construct from Istream
template<class ReactionThermo>
Reaction<ReactionThermo>::Reaction
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    Istream& is
)
:
    ReactionThermo(*thermoDatabase[species[0]]),
    species_(species)
{
    setLRhs(is);
    setThermo(thermoDatabase);
}
    

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class ReactionThermo>
autoPtr<Reaction<ReactionThermo> > Reaction<ReactionThermo>::New
(
    const speciesTable& species,
    const HashPtrTable<ReactionThermo>& thermoDatabase,
    Istream& is
)
{
    if (is.eof())
    {
        FatalIOErrorIn
        (
            "Reaction<ReactionThermo>::New(const speciesTable& species,"
            " const HashPtrTable<ReactionThermo>& thermoDatabase, Istream&)",
            is
        )   << "Reaction type not specified" << endl << endl
            << "Valid Reaction types are :" << endl
            << IstreamConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    word reactionTypeName(is);

    typename IstreamConstructorTable::iterator cstrIter
        = IstreamConstructorTablePtr_->find(reactionTypeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "Reaction<ReactionThermo>::New(const speciesTable& species,"
            " const HashPtrTable<ReactionThermo>& thermoDatabase, Istream&)",
            is
        )   << "Unknown reaction type " << reactionTypeName << endl << endl
            << "Valid reaction types are :" << endl
            << IstreamConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return autoPtr<Reaction<ReactionThermo> >
    (
        cstrIter()(species, thermoDatabase, is)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ReactionThermo>
void Reaction<ReactionThermo>::write(Ostream& os) const
{
    os << type() << nl << "    ";

    forAll(lhs_, i)
    {
        const typename Reaction<ReactionThermo>::specieCoeffs& sc = lhs_[i];

        if (sc.stoichCoeff != 1)
        {
            os << sc.stoichCoeff;
        }

        os << species_[sc.index];

        if (sc.exponent != sc.stoichCoeff)
        {
            os << '^' << sc.exponent;
        }

        if (i < lhs_.size() - 1)
        {
            os << " + ";
        }
    }

    os << " = ";

    forAll(rhs_, i)
    {
        const typename Reaction<ReactionThermo>::specieCoeffs& sc = rhs_[i];

        if (sc.stoichCoeff != 1)
        {
            os << sc.stoichCoeff;
        }

        os << species_[sc.index];

        if (sc.exponent != sc.stoichCoeff)
        {
            os << '^' << sc.exponent;
        }

        if (i < rhs_.size() - 1)
        {
            os << " + ";
        }
    }

    os  << endl << "   ";
}


template<class ReactionThermo>
scalar Reaction<ReactionThermo>::kf
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return 0.0;
}


template<class ReactionThermo>
scalar Reaction<ReactionThermo>::kr
(
    const scalar kfwd,
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return 0.0;
}

template<class ReactionThermo>
scalar Reaction<ReactionThermo>::kr
(
    const scalar T,
    const scalar p,
    const scalarField& c
) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
