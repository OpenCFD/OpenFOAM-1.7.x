/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::referredMolecule::id() const
{
    return id_;
}


inline const Foam::vector& Foam::referredMolecule::position() const
{
    return position_;
}


inline const Foam::List<Foam::vector>&
Foam::referredMolecule::sitePositions() const
{
    return sitePositions_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==
(
    const Foam::referredMolecule& a,
    const Foam::referredMolecule& b
)
{
    return
    (
        a.id() == b.id()
     && a.position() == b.position()
    );
}


inline bool Foam::operator!=
(
    const Foam::referredMolecule& a,
    const Foam::referredMolecule& b
)
{
    return !(a == b);
}


// ************************************************************************* //
