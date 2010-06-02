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

#include "boundBox.H"
#include "PstreamReduceOps.H"
#include "tmp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::scalar Foam::boundBox::great(VGREAT);

const Foam::boundBox Foam::boundBox::greatBox
(
    point(-VGREAT, -VGREAT, -VGREAT),
    point(VGREAT, VGREAT, VGREAT)
);


const Foam::boundBox Foam::boundBox::invertedBox
(
    point(VGREAT, VGREAT, VGREAT),
    point(-VGREAT, -VGREAT, -VGREAT)
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boundBox::calculate(const pointField& points, const bool doReduce)
{
    if (points.empty())
    {
        min_ = point::zero;
        max_ = point::zero;

        if (doReduce && Pstream::parRun())
        {
            // Use values that get overwritten by reduce minOp, maxOp below
            min_ = point(VGREAT, VGREAT, VGREAT);
            max_ = point(-VGREAT, -VGREAT, -VGREAT);
        }
    }
    else
    {
        min_ = points[0];
        max_ = points[0];

        for (label i = 1; i < points.size(); i++)
        {
            min_ = ::Foam::min(min_, points[i]);
            max_ = ::Foam::max(max_, points[i]);
        }
    }

    // Reduce parallel information
    if (doReduce)
    {
        reduce(min_, minOp<point>());
        reduce(max_, maxOp<point>());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundBox::boundBox(const pointField& points, const bool doReduce)
:
    min_(point::zero),
    max_(point::zero)
{
    calculate(points, doReduce);
}


Foam::boundBox::boundBox(const tmp<pointField>& points, const bool doReduce)
:
    min_(point::zero),
    max_(point::zero)
{
    calculate(points(), doReduce);
    points.clear();
}


Foam::boundBox::boundBox(Istream& is)
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const boundBox& bb)
{
    if (os.format() == IOstream::ASCII)
    {
        os << bb.min_ << token::SPACE << bb.max_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&bb.min_),
            sizeof(boundBox)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const boundBox&)");
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, boundBox& bb)
{
    if (is.format() == IOstream::ASCII)
    {
        return is >> bb.min_ >> bb.max_;
    }
    else
    {
        is.read
        (
            reinterpret_cast<char*>(&bb.min_),
            sizeof(boundBox)
        );
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, boundBox&)");
    return is;
}

// ************************************************************************* //
