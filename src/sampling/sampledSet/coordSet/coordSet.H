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

Class
    Foam::coordSet

Description
    Holds list of sampling positions

SourceFiles
    coordSet.C

\*---------------------------------------------------------------------------*/

#ifndef coordSet_H
#define coordSet_H

#include "pointField.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class coordSet Declaration
\*---------------------------------------------------------------------------*/

class coordSet
:
    public pointField
{

protected:

        //- Name
        const word name_;

        //- Axis write type
        const word axis_;

        //- Reference point for "distance" write specifier.
        point refPoint_;

public:

    // Constructors

        //- Construct from components
        coordSet
        (
            const word& name,
            const word& axis
        );


        //- Construct from components
        coordSet
        (
            const word& name,
            const word& axis,
            const List<point>& points,
            const point& refPoint = point::zero
        );

        //- Construct from components
        coordSet
        (
            const word& name,
            const word& axis,
            const scalarField& points,
            const scalar refPoint = 0.0
        );


    // Member functions

        const word& name() const
        {
            return name_;
        }

        const word& axis() const
        {
            return axis_;
        }

        const point& refPoint() const
        {
            return refPoint_;
        }

        //- Is axis specification a vector
        bool hasVectorAxis() const;

        //- Get coordinate of point according to axis specification.
        //  If axis="distance" can be: -distance to starting point (e.g.
        //  uniformSet) or -distance to first sampling point
        //  (e.g. cloudSet)
        scalar scalarCoord
        (
            const label index
        ) const;

        //- Get point according to axis="full" specification
        vector vectorCoord
        (
            const label index
        ) const;

        Ostream& write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
