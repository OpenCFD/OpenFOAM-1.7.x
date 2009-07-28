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

\*---------------------------------------------------------------------------*/

#include "icoPolynomial.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<int PolySize>
inline Foam::icoPolynomial<PolySize>::icoPolynomial
(
    const specie& sp,
    const Polynomial<PolySize>& rhoPoly
)
:
    specie(sp),
    rhoPolynomial_(rhoPoly)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int PolySize>
inline Foam::icoPolynomial<PolySize>::icoPolynomial
(
    const word& name,
    const icoPolynomial<PolySize>& ip
)
:
    specie(name, ip),
    rhoPolynomial_(ip.rhoPolynomial_)
{}


template<int PolySize>
inline Foam::autoPtr<Foam::icoPolynomial<PolySize> >
Foam::icoPolynomial<PolySize>::clone() const
{
    return autoPtr<icoPolynomial<PolySize> >
    (
        new icoPolynomial<PolySize>(*this)
    );
}


template<int PolySize>
inline Foam::autoPtr<Foam::icoPolynomial<PolySize> >
Foam::icoPolynomial<PolySize>::New(Istream& is)
{
    return autoPtr<icoPolynomial<PolySize> >(new icoPolynomial<PolySize>(is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int PolySize>
inline Foam::scalar Foam::icoPolynomial<PolySize>::rho(scalar, scalar T) const
{
    return rhoPolynomial_.evaluate(T);
}


template<int PolySize>
inline Foam::scalar Foam::icoPolynomial<PolySize>::psi(scalar, scalar) const
{
    return 0.0;
}


template<int PolySize>
inline Foam::scalar Foam::icoPolynomial<PolySize>::Z(scalar, scalar) const
{
    return 0.0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<int PolySize>
inline void Foam::icoPolynomial<PolySize>::operator+=
(
    const icoPolynomial<PolySize>& ip
)
{
    scalar molr1 = this->nMoles();

    specie::operator+=(ip);

    molr1 /= this->nMoles();
    scalar molr2 = ip.nMoles()/this->nMoles();

    rhoPolynomial_ = molr1*rhoPolynomial_ + molr2*ip.rhoPolynomial_;
}


template<int PolySize>
inline void Foam::icoPolynomial<PolySize>::operator-=
(
    const icoPolynomial<PolySize>& ip
)
{
    scalar molr1 = this->nMoles();
    specie::operator-=(ip);

    molr1 /= this->nMoles();
    scalar molr2 = ip.nMoles()/this->nMoles();

    rhoPolynomial_ = molr1*rhoPolynomial_ - molr2*ip.rhoPolynomial_;
}


template<int PolySize>
inline void Foam::icoPolynomial<PolySize>::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<int PolySize>
Foam::icoPolynomial<PolySize> Foam::operator+
(
    const icoPolynomial<PolySize>& ip1,
    const icoPolynomial<PolySize>& ip2
)
{
    scalar mol1 = ip1.nMoles();
    scalar mol2 = ip2.nMoles();
    scalar nMoles = mol1 + mol2;

    return icoPolynomial<PolySize>
    (
        static_cast<const specie&>(ip1)
      + static_cast<const specie&>(ip2),
        (mol1/nMoles)*ip1.rhoPolynomial_ + (mol2/nMoles)*ip2.rhoPolynomial_
    );
}


template<int PolySize>
Foam::icoPolynomial<PolySize> Foam::operator-
(
    const icoPolynomial<PolySize>& ip1,
    const icoPolynomial<PolySize>& ip2
)
{
    scalar mol1 = ip1.nMoles();
    scalar mol2 = ip2.nMoles();
    scalar nMoles = mol1 + mol2;

    return icoPolynomial<PolySize>
    (
        static_cast<const specie&>(ip1)
      - static_cast<const specie&>(ip2),
        (mol1/nMoles)*ip1.rhoPolynomial_ - (mol2/nMoles)*ip2.rhoPolynomial_
    );
}


template<int PolySize>
Foam::icoPolynomial<PolySize> Foam::operator*
(
    const scalar s,
    const icoPolynomial<PolySize>& ip
)
{
    return icoPolynomial<PolySize>
    (
        s*static_cast<const specie&>(ip),
        ip.rhoPolynomial_
    );
}


template<int PolySize>
Foam::icoPolynomial<PolySize> Foam::operator==
(
    const icoPolynomial<PolySize>& ip1,
    const icoPolynomial<PolySize>& ip2
)
{
    return ip2 - ip1;
}


// ************************************************************************* //
