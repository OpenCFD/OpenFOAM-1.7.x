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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class equationOfState>
inline Foam::eConstThermo<equationOfState>::eConstThermo
(
    const equationOfState& st,
    const scalar cv,
    const scalar hf
)
:
    equationOfState(st),
    Cv_(cv),
    Hf_(hf)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::eConstThermo<equationOfState>::eConstThermo
(
    const word& name,
    const eConstThermo<equationOfState>& ct
)
:
    equationOfState(name, ct),
    Cv_(ct.Cv_),
    Hf_(ct.Hf_)
{}


template<class equationOfState>
inline Foam::autoPtr<Foam::eConstThermo<equationOfState> >
Foam::eConstThermo<equationOfState>::clone() const
{
    return autoPtr<eConstThermo<equationOfState> >
    (
        new eConstThermo<equationOfState>(*this)
    );
}


template<class equationOfState>
inline Foam::autoPtr<Foam::eConstThermo<equationOfState> >
Foam::eConstThermo<equationOfState>::New(Istream& is)
{
    return autoPtr<eConstThermo<equationOfState> >
    (
        new eConstThermo<equationOfState>(is)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::scalar Foam::eConstThermo<equationOfState>::cp
(
    const scalar
) const
{
    return Cv_*this->W() + specie::RR;
}


template<class equationOfState>
inline Foam::scalar Foam::eConstThermo<equationOfState>::h
(
    const scalar T
) const
{
    return cp(T)*T + Hf_*this->W();
}


template<class equationOfState>
inline Foam::scalar Foam::eConstThermo<equationOfState>::hs
(
    const scalar T
) const
{
    return cp(T)*T;
}


template<class equationOfState>
inline Foam::scalar Foam::eConstThermo<equationOfState>::hc() const
{
    return Hf_*this->W();
}


template<class equationOfState>
inline Foam::scalar Foam::eConstThermo<equationOfState>::s
(
    const scalar T
) const
{
    notImplemented("scalar eConstThermo::s(const scalar T) const");
    return T;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline void Foam::eConstThermo<equationOfState>::operator+=
(
    const eConstThermo<equationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    equationOfState::operator+=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Cv_ = molr1*Cv_ + molr2*ct.Cv_;
    Hf_ = molr1*Hf_ + molr2*ct.Hf_;
}


template<class equationOfState>
inline void Foam::eConstThermo<equationOfState>::operator-=
(
    const eConstThermo<equationOfState>& ct
)
{
    scalar molr1 = this->nMoles();

    equationOfState::operator-=(ct);

    molr1 /= this->nMoles();
    scalar molr2 = ct.nMoles()/this->nMoles();

    Cv_ = molr1*Cv_ - molr2*ct.Cv_;
    Hf_ = molr1*Hf_ - molr2*ct.Hf_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class equationOfState>
inline Foam::eConstThermo<equationOfState> Foam::operator+
(
    const eConstThermo<equationOfState>& ct1,
    const eConstThermo<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
      + static_cast<const equationOfState&>(ct2)
    );

    return eConstThermo<equationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.Cv_
      + ct2.nMoles()/eofs.nMoles()*ct2.Cv_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      + ct2.nMoles()/eofs.nMoles()*ct2.Hf_
    );
}


template<class equationOfState>
inline Foam::eConstThermo<equationOfState> Foam::operator-
(
    const eConstThermo<equationOfState>& ct1,
    const eConstThermo<equationOfState>& ct2
)
{
    equationOfState eofs
    (
        static_cast<const equationOfState&>(ct1)
      - static_cast<const equationOfState&>(ct2)
    );

    return eConstThermo<equationOfState>
    (
        eofs,
        ct1.nMoles()/eofs.nMoles()*ct1.Cv_
      - ct2.nMoles()/eofs.nMoles()*ct2.Cv_,
        ct1.nMoles()/eofs.nMoles()*ct1.Hf_
      - ct2.nMoles()/eofs.nMoles()*ct2.Hf_
    );
}


template<class equationOfState>
inline Foam::eConstThermo<equationOfState> Foam::operator*
(
    const scalar s,
    const eConstThermo<equationOfState>& ct
)
{
    return eConstThermo<equationOfState>
    (
        s*static_cast<const equationOfState&>(ct),
        ct.Cv_,
        ct.Hf_
    );
}


template<class equationOfState>
inline Foam::eConstThermo<equationOfState> Foam::operator==
(
    const eConstThermo<equationOfState>& ct1,
    const eConstThermo<equationOfState>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
