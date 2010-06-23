/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "basicSource.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::basicSource::name() const
{
    return name_;
}


inline const Foam::fvMesh& Foam::basicSource::mesh() const
{
    return mesh_;
}

inline const Foam::dictionary& Foam::basicSource::dictCoeffs() const
{
    return dict_;
}

inline bool Foam::basicSource::active() const
{
    return active_;
}


inline Foam::scalar Foam::basicSource::timeStart() const
{
    return timeStart_;
}


inline Foam::scalar Foam::basicSource::duration() const
{
    return duration_;
}


inline Foam::scalar Foam::basicSource::timeEnd() const
{
    return timeStart_ + duration_;
}


inline const Foam::basicSource::selectionModeType&
Foam::basicSource::selectionMode() const
{
    return selectionMode_;
}


inline const Foam::word& Foam::basicSource::cellSetName() const
{
    return cellSetName_;
}


inline Foam::scalar Foam::basicSource::V() const
{
    return V_;
}


inline const Foam::labelList& Foam::basicSource::cells() const
{
    return cells_;
}


inline Foam::word& Foam::basicSource::name()
{
    return name_;
}


inline bool& Foam::basicSource::active()
{
    return active_;
}


inline Foam::scalar& Foam::basicSource::timeStart()
{
    return timeStart_;
}


inline Foam::scalar& Foam::basicSource::duration()
{
    return duration_;
}


inline Foam::word& Foam::basicSource::cellSetName()
{
    return cellSetName_;
}


inline Foam::basicSource::selectionModeType&
Foam::basicSource::selectionMode()
{
    return selectionMode_;
}


inline Foam::scalar& Foam::basicSource::V()
{
    return V_;
}


inline Foam::labelList& Foam::basicSource::cells()
{
    return cells_;
}


// ************************************************************************* //
