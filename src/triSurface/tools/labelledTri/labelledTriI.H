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

Description

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
inline labelledTri::labelledTri()
:
    region_(-1)
{}


//- Construct from components
inline labelledTri::labelledTri
(
    const label A,
    const label B,
    const label C,
    const label region
)
:
    triFace(A, B, C),
    region_(region)
{}


inline labelledTri::labelledTri(Istream& is)
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline label labelledTri::region() const
{
    return region_;
}

inline label& labelledTri::region()
{
    return region_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Istream& operator>>(Istream& is, labelledTri& t)
{
    if (is.format() == IOstream::ASCII)
    {
        // Read beginning of labelledTri point pair
        is.readBegin("labelledTri");

        is  >> static_cast<triFace&>(t) >> t.region_;

        // Read end of labelledTri point pair
        is.readEnd("labelledTri");
    }
    else
    {
        is.read(reinterpret_cast<char*>(&t), sizeof(labelledTri));
    }

    // Check state of Ostream
    is.check("Istream& operator>>(Istream&, labelledTri&)");

    return is;
}


inline Ostream& operator<<(Ostream& os, const labelledTri& t)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << token::BEGIN_LIST
            << static_cast<const triFace&>(t) << token::SPACE << t.region_
            << token::END_LIST;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&t),
            sizeof(labelledTri)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const labelledTri&)");


    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
