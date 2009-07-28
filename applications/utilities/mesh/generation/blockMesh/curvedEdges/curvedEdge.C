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
    library functions that will define a curvedEdge in space
    parameterised for 0<lambda<1 from the beginning
    point to the end point.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "curvedEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(curvedEdge, 0);
defineRunTimeSelectionTable(curvedEdge, Istream);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
curvedEdge::curvedEdge
(
    const pointField& points,
    const label start,
    const label end
)
:
    points_(points),
    start_(start),
    end_(end)
{}


// Construct from Istream
curvedEdge::curvedEdge(const pointField& points, Istream& is)
:
    points_(points),
    start_(readLabel(is)),
    end_(readLabel(is))
{}


// Copy construct
curvedEdge::curvedEdge(const curvedEdge& c)
:
    points_(c.points_),
    start_(c.start_),
    end_(c.end_)
{}


//- Clone function
autoPtr<curvedEdge> curvedEdge::clone() const
{
    notImplemented("curvedEdge::clone() const");
    return autoPtr<curvedEdge>(NULL);
}


//- New function which constructs and returns pointer to a curvedEdge
autoPtr<curvedEdge> curvedEdge::New(const pointField& points, Istream& is)
{
    if (debug)
    {
        Info<< "curvedEdge::New(const pointField&, Istream&) : "
            << "constructing curvedEdge"
            << endl;
    }

    word curvedEdgeType(is);

    IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_
            ->find(curvedEdgeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorIn("curvedEdge::New(const pointField&, Istream&)")
            << "Unknown curvedEdge type " << curvedEdgeType << endl << endl
            << "Valid curvedEdge types are" << endl
            << IstreamConstructorTablePtr_->toc()
            << abort(FatalError);
    }

    return autoPtr<curvedEdge>(cstrIter()(points, is));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return the complete knotList by adding the start and end points to the
//  given list
pointField knotlist
(
    const pointField& points,
    const label start,
    const label end,
    const pointField& otherknots
)
{
    label listsize(otherknots.size() + 2);
    pointField tmp(listsize);

    tmp[0] = points[start];

    for (register label i=1; i<listsize-1; i++)
    {
        tmp[i] = otherknots[i-1];
    }

    tmp[listsize-1] = points[end];

    return tmp;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void curvedEdge::operator=(const curvedEdge&)
{
    notImplemented("void curvedEdge::operator=(const curvedEdge&)");
}


Ostream& operator<<(Ostream& os, const curvedEdge& p)
{
    os << p.start_ << tab << p.end_ << endl;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
