/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
inline LangmuirHinshelwoodReactionRate::LangmuirHinshelwoodReactionRate
(
    const scalar A[],
    const scalar Ta[],
    const label co,
    const label c3h6,
    const label no
)
:
    co_(co),
    c3h6_(c3h6),
    no_(no)
{
    for (int i=0; i<n_; i++)
    {
        A_[i] = A[i];
        Ta_[i] = Ta[i];
    }
}


//- Construct from Istream
inline LangmuirHinshelwoodReactionRate::LangmuirHinshelwoodReactionRate
(
    const speciesTable& st,
    Istream& is
)
:
    co_(st["CO"]),
    c3h6_(st["C3H6"]),
    no_(st["NO"])
{
    is.readBegin("LangmuirHinshelwoodReactionRate(Istream&)");

    for (int i=0; i<n_; i++)
    {
        is >> A_[i] >> Ta_[i];
    }

    is.readEnd("LangmuirHinshelwoodReactionRate(Istream&)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar LangmuirHinshelwoodReactionRate::operator()
(
    const scalar T,
    const scalar,
    const scalarField& c
) const
{
    return A_[0]*exp(-Ta_[0]/T)/
    (
        T
       *sqr(1 + A_[1]*exp(-Ta_[1]/T)*c[co_] + A_[2]*exp(-Ta_[2]/T)*c[c3h6_])
       *(1 + A_[3]*exp(-Ta_[3]/T)*sqr(c[co_])*sqr(c[c3h6_]))
       *(1 + A_[4]*exp(-Ta_[4]/T)*pow(c[no_], 0.7))
    );
}


inline Ostream& operator<<
(
    Ostream& os,
    const LangmuirHinshelwoodReactionRate& lhrr
)
{
    os  << token::BEGIN_LIST;

    for (int i=0; i<LangmuirHinshelwoodReactionRate::n_; i++)
    {
        os  << token::SPACE << lhrr.A_[i] << token::SPACE << lhrr.Ta_[i];
    }

    os << token::END_LIST;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
