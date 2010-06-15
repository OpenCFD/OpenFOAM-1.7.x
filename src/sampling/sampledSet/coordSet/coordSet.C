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

#include "coordSet.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
Foam::coordSet::coordSet
(
    const word& name,
    const word& axis
)
:
    pointField(0),
    name_(name),
    axis_(axis),
    refPoint_(vector::zero)
{}


//- Construct from components
Foam::coordSet::coordSet
(
    const word& name,
    const word& axis,
    const List<point>& points,
    const point& refPoint
)
:
    pointField(points),
    name_(name),
    axis_(axis),
    refPoint_(refPoint)
{}


//- Construct from components
Foam::coordSet::coordSet
(
    const word& name,
    const word& axis,
    const scalarField& points,
    const scalar refPoint
)
:
    pointField(points.size(), point::zero),
    name_(name),
    axis_(axis),
    refPoint_(point::zero)
{
    if (axis_ == "x" || axis_ == "distance")
    {
        refPoint_.x() = refPoint;
        replace(point::X, points);
    }
    else if (axis_ == "y")
    {
        replace(point::Y, points);
    }
    else if (axis_ == "z")
    {
        replace(point::Z, points);
    }
    else
    {
        FatalErrorIn
        (
            "coordSet::coordSet(const word& name,"
            "const word& axis, const List<scalar>& points,"
            "const scalar refPoint)"
        )   << "Illegal axis specification " << axis_
            << " for sampling line " << name_
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::coordSet::hasVectorAxis() const
{
    return axis_ == "xyz";
}


Foam::scalar Foam::coordSet::scalarCoord
(
    const label index
)   const
{
    const point& p = operator[](index);

    if (axis_ == "x")
    {
        return p.x();
    }
    else if (axis_ == "y")
    {
        return p.y();
    }
    else if (axis_ == "z")
    {
        return p.z();
    }
    else if (axis_ == "distance")
    {
        // Use distance to reference point
        return mag(p - refPoint_);
    }
    else
    {
        FatalErrorIn
        (
            "coordSet::scalarCoord(const label)"
        )   << "Illegal axis specification " << axis_
            << " for sampling line " << name_
            << exit(FatalError);

        return 0;
    }
}


Foam::point Foam::coordSet::vectorCoord(const label index) const
{
    const point& p = operator[](index);

    return p;
}


Foam::Ostream& Foam::coordSet::write(Ostream& os) const
{
    os  << "name:" << name_ << " axis:" << axis_ << " reference:" << refPoint_
        << endl
        << endl << "\t(coord)"
        << endl;

    forAll(*this, sampleI)
    {
        os  << '\t' << operator[](sampleI) << endl;
    }

    return os;
}


// ************************************************************************* //
