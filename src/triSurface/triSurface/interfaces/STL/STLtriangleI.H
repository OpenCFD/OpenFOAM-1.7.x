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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
inline STLtriangle::STLtriangle()
{}


//- Construct from components
inline STLtriangle::STLtriangle
(
    const STLpoint& normal,
    const STLpoint& a,
    const STLpoint& b,
    const STLpoint& c,
    unsigned short region
)
:
    normal_(normal),
    a_(a),
    b_(b),
    c_(c),
    region_(region)
{}


//- Construct from istream (binary)
inline STLtriangle::STLtriangle(istream& is)
{
    read(is);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const STLpoint& STLtriangle::a() const
{
    return a_;
}

inline const STLpoint& STLtriangle::b() const
{
    return b_;
}

inline const STLpoint& STLtriangle::c() const
{
    return c_;
}

inline unsigned short STLtriangle::region() const
{
    return region_;
}


inline void STLtriangle::read(istream& is)
{
    is.read(reinterpret_cast<char*>(this), 4*sizeof(STLpoint));
    is.read(reinterpret_cast<char*>(&region_), 2);
}


inline void STLtriangle::write(ostream& os)
{
    os.write(reinterpret_cast<char*>(this), 4*sizeof(STLpoint));
    os.write(reinterpret_cast<char*>(&region_), 2);
}

// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Ostream& operator<<(Ostream& os, const STLtriangle& stlt)
{
    os  << stlt.normal_ << token::SPACE
        << stlt.a_ << token::SPACE
        << stlt.b_ << token::SPACE
        << stlt.c_ << token::SPACE
        << stlt.region_;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
