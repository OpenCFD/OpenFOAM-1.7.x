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

#include "Polynomial.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<int PolySize>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Polynomial<PolySize>& poly
)
{
    os  << static_cast
            <VectorSpace<Polynomial<PolySize>, scalar, PolySize> >(poly);

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const Polynomial<PolySize>&)"
    );

    return os;
}


// ************************************************************************* //
