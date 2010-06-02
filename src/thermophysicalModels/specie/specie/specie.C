/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
#include "IOstreams.H"
#include "dimensionedConstants.H"

/* * * * * * * * * * * * * public constants  * * * * * * * * * * * * */

//- Universal gas constant (default in [J/(kmol K)])
const Foam::scalar Foam::specie::RR = dimensionedConstant("R", 8314.51);

//- Standard pressure (default in [Pa])
const Foam::scalar Foam::specie::Pstd = dimensionedConstant("Pstd", 1.0e5);

//- Standard temperature (default in [K])
const Foam::scalar Foam::specie::Tstd = dimensionedConstant("Tstd", 298.15);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specie::specie(Istream& is)
:
    name_(is),
    nMoles_(readScalar(is)),
    molWeight_(readScalar(is))
{
    is.check("specie::specie(Istream& is)");
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const specie& st)
{
    os  << st.name_ << tab
        << st.nMoles_ << tab
        << st.molWeight_;

    os.check("Ostream& operator<<(Ostream& os, const specie& st)");
    return os;
}


// ************************************************************************* //
