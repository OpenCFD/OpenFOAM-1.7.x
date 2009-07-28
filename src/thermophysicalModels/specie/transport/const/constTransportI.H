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
template<class thermo>
inline constTransport<thermo>::constTransport
(
    const thermo& t,
    const scalar mu,
    const scalar Pr
)
:
    thermo(t),
    Mu(mu),
    rPr(1.0/Pr)
{}


// Construct as named copy
template<class thermo>
inline constTransport<thermo>::constTransport
(
    const word& name,
    const constTransport& ct
)
:
    thermo(name, ct),
    Mu(ct.Mu),
    rPr(ct.rPr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Dynamic viscosity [kg/ms]
template<class thermo>
inline scalar constTransport<thermo>::mu(const scalar) const
{
    return Mu;
}


// Thermal conductivity [W/mK]
template<class thermo>
inline scalar constTransport<thermo>::kappa(const scalar T) const
{
    return this->Cp(T)*mu(T)*rPr;
}


// Thermal diffusivity for enthalpy [kg/ms]
template<class thermo>
inline scalar constTransport<thermo>::alpha(const scalar T) const
{
    scalar Cp_ = this->Cp(T);

    scalar deltaT = T - specie::Tstd;
    scalar CpBar =
        (deltaT*(this->H(T) - this->H(specie::Tstd)) + Cp_)/(sqr(deltaT) + 1);

    return Cp_*mu(T)*rPr/CpBar;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class thermo>
inline constTransport<thermo>& constTransport<thermo>::operator=
(
    const constTransport<thermo>& ct
)
{
    thermo::operator=(ct);

    Mu = ct.Mu;
    rPr = ct.rPr;

    return *this;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class thermo>
inline constTransport<thermo> operator+
(
    const constTransport<thermo>& ct1,
    const constTransport<thermo>& ct2
)
{
    thermo t
    (
        static_cast<const thermo&>(ct1) + static_cast<const thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return constTransport<thermo>
    (
        t,
        molr1*ct1.Mu + molr2*ct2.Mu,
        molr1*ct1.rPr + molr2*ct2.rPr
    );
}


template<class thermo>
inline constTransport<thermo> operator-
(
    const constTransport<thermo>& ct1,
    const constTransport<thermo>& ct2
)
{
    thermo t
    (
        static_cast<const thermo&>(ct1) - static_cast<const thermo&>(ct2)
    );

    scalar molr1 = ct1.nMoles()/t.nMoles();
    scalar molr2 = ct2.nMoles()/t.nMoles();

    return constTransport<thermo>
    (
        t,
        molr1*ct1.Mu - molr2*ct2.Mu,
        molr1*ct1.rPr - molr2*ct2.rPr
    );
}


template<class thermo>
inline constTransport<thermo> operator*
(
    const scalar s,
    const constTransport<thermo>& ct
)
{
    return constTransport<thermo>
    (
        s*static_cast<const thermo&>(ct),
        ct.Mu,
        ct.rPr
    );
}


template<class thermo>
inline constTransport<thermo> operator==
(
    const constTransport<thermo>& ct1,
    const constTransport<thermo>& ct2
)
{
    return ct2 - ct1;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
