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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const char* const pTraits<Scalar>::typeName = "scalar";
const Scalar pTraits<Scalar>::zero = 0.0;
const Scalar pTraits<Scalar>::one = 1.0;
const Scalar pTraits<Scalar>::min = -ScalarVGREAT;
const Scalar pTraits<Scalar>::max = ScalarVGREAT;

const char* pTraits<Scalar>::componentNames[] = { "x" };

pTraits<Scalar>::pTraits(Istream& is)
{
    is >> p_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

word name(const Scalar val)
{
    std::ostringstream buf;
    buf << val;
    return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Scalar readScalar(Istream& is)
{
    Scalar rs;
    is >> rs;

    return rs;
}


Istream& operator>>(Istream& is, Scalar& s)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isNumber())
    {
        s = t.number();
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, Scalar&)", is)
            << "wrong token type - expected Scalar found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, Scalar&)");

    return is;
}


Ostream& operator<<(Ostream& os, const Scalar s)
{
    os.write(s);
    os.check("Ostream& operator<<(Ostream&, const Scalar&)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
