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

\*---------------------------------------------------------------------------*/

#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::polynomialTransport<Thermo, PolySize>::polynomialTransport
(
    const polynomialTransport& pt
)
:
    Thermo(pt),
    muPolynomial_(pt.muPolynomial_),
    kappaPolynomial_(pt.kappaPolynomial_)
{}


template<class Thermo, int PolySize>
inline Foam::polynomialTransport<Thermo, PolySize>::polynomialTransport
(
    const Thermo& t,
    const Polynomial<PolySize>& muPoly,
    const Polynomial<PolySize>& kappaPoly
)
:
    Thermo(t),
    muPolynomial_(muPoly),
    kappaPolynomial_(kappaPoly)
{}


template<class Thermo, int PolySize>
inline Foam::polynomialTransport<Thermo, PolySize>::polynomialTransport
(
    const word& name,
    const polynomialTransport& pt
)
:
    Thermo(name, pt),
    muPolynomial_(pt.muPolynomial_),
    kappaPolynomial_(pt.kappaPolynomial_)
{}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::polynomialTransport<Thermo, PolySize> >
Foam::polynomialTransport<Thermo, PolySize>::clone() const
{
    return autoPtr<polynomialTransport<Thermo, PolySize> >
    (
        new polynomialTransport<Thermo, PolySize>(*this)
    );
}


template<class Thermo, int PolySize>
inline Foam::autoPtr<Foam::polynomialTransport<Thermo, PolySize> >
Foam::polynomialTransport<Thermo, PolySize>::New(Istream& is)
{
    return autoPtr<polynomialTransport<Thermo, PolySize> >
    (
        new polynomialTransport<Thermo, PolySize>(is)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::scalar Foam::polynomialTransport<Thermo, PolySize>::mu
(
    const scalar T
) const
{
    return muPolynomial_.evaluate(T)/this->W();
}


template<class Thermo, int PolySize>
inline Foam::scalar Foam::polynomialTransport<Thermo, PolySize>::kappa
(
    const scalar T
) const
{
    return kappaPolynomial_.evaluate(T)/this->W();
}


template<class Thermo, int PolySize>
inline Foam::scalar Foam::polynomialTransport<Thermo, PolySize>::alpha
(
    const scalar T
) const
{
    scalar deltaT = T - specie::Tstd;
    scalar CpBar =
        (deltaT*(this->H(T) - this->H(specie::Tstd)) + this->Cp(T))
       /(sqr(deltaT) + 1);

    return kappa(T)/CpBar;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::polynomialTransport<Thermo, PolySize>&
Foam::polynomialTransport<Thermo, PolySize>::operator=
(
    const polynomialTransport<Thermo, PolySize>& pt
)
{
    Thermo::operator=(pt);

    muPolynomial_ = pt.muPolynomial_;
    kappaPolynomial_ = pt.kappaPolynomial_;

    return *this;
}


template<class Thermo, int PolySize>
inline void Foam::polynomialTransport<Thermo, PolySize>::operator+=
(
    const polynomialTransport<Thermo, PolySize>& pt
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator+=(pt);

    molr1 /= this->nMoles();
    scalar molr2 = pt.nMoles()/this->nMoles();

    muPolynomial_ = molr1*muPolynomial_ + molr2*pt.muPolynomial_;
    kappaPolynomial_ = molr1*kappaPolynomial_ + molr2*pt.kappaPolynomial_;
}


template<class Thermo, int PolySize>
inline void Foam::polynomialTransport<Thermo, PolySize>::operator-=
(
    const polynomialTransport<Thermo, PolySize>& pt
)
{
    scalar molr1 = this->nMoles();

    Thermo::operator-=(pt);

    molr1 /= this->nMoles();
    scalar molr2 = pt.nMoles()/this->nMoles();

    muPolynomial_ = molr1*muPolynomial_ - molr2*pt.muPolynomial_;
    kappaPolynomial_ = molr1*kappaPolynomial_ - molr2*pt.kappaPolynomial_;
}


template<class Thermo, int PolySize>
inline void Foam::polynomialTransport<Thermo, PolySize>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
inline Foam::polynomialTransport<Thermo, PolySize> Foam::operator+
(
    const polynomialTransport<Thermo, PolySize>& pt1,
    const polynomialTransport<Thermo, PolySize>& pt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(pt1) + static_cast<const Thermo&>(pt2)
    );

    scalar molr1 = pt1.nMoles()/t.nMoles();
    scalar molr2 = pt2.nMoles()/t.nMoles();

    return polynomialTransport<Thermo, PolySize>
    (
        t,
        molr1*pt1.muPolynomial_ + molr2*pt2.muPolynomial_,
        molr1*pt1.kappaPolynomial_ + molr2*pt2.kappaPolynomial_
    );
}


template<class Thermo, int PolySize>
inline Foam::polynomialTransport<Thermo, PolySize> Foam::operator-
(
    const polynomialTransport<Thermo, PolySize>& pt1,
    const polynomialTransport<Thermo, PolySize>& pt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(pt1) - static_cast<const Thermo&>(pt2)
    );

    scalar molr1 = pt1.nMoles()/t.nMoles();
    scalar molr2 = pt2.nMoles()/t.nMoles();

    return polynomialTransport<Thermo, PolySize>
    (
        t,
        molr1*pt1.muPolynomial_ - molr2*pt2.muPolynomial_,
        molr1*pt1.kappaPolynomial_ - molr2*pt2.kappaPolynomial_
    );
}


template<class Thermo, int PolySize>
inline Foam::polynomialTransport<Thermo, PolySize> Foam::operator*
(
    const scalar s,
    const polynomialTransport<Thermo, PolySize>& pt
)
{
    return polynomialTransport<Thermo, PolySize>
    (
        s*static_cast<const Thermo&>(pt),
        pt.muPolynomial_,
        pt.kappaPolynomial_
    );
}


template<class Thermo, int PolySize>
inline Foam::polynomialTransport<Thermo, PolySize> Foam::operator==
(
    const polynomialTransport<Thermo, PolySize>& pt1,
    const polynomialTransport<Thermo, PolySize>& pt2
)
{
    return pt2 - pt1;
}


// ************************************************************************* //
