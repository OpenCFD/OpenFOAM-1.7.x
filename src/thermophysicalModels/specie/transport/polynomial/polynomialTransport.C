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

#include "polynomialTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::polynomialTransport<Thermo, PolySize>::polynomialTransport(Istream& is)
:
    Thermo(is),
    muPolynomial_("muPolynomial", is),
    kappaPolynomial_("kappaPolynomial", is)
{
    muPolynomial_ *= this->W();
    kappaPolynomial_ *= this->W();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo, int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const polynomialTransport<Thermo, PolySize>& pt
)
{
    os  << static_cast<const Thermo&>(pt) << tab
        << "muPolynomial" << tab << pt.muPolynomial_/pt.W() << tab
        << "kappaPolynomial" << tab << pt.kappaPolynomial_/pt.W();

    os.check
    (
        "Ostream& operator<<"
        "("
            "Ostream&, "
            "const polynomialTransport<Thermo, PolySize>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
