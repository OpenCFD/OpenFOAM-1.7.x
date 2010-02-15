/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "kinematicParcelInjectionData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::point& Foam::kinematicParcelInjectionData::x() const
{
    return x_;
}


inline const Foam::vector& Foam::kinematicParcelInjectionData::U() const
{
    return U_;
}


inline Foam::scalar Foam::kinematicParcelInjectionData::d() const
{
    return d_;
}


inline Foam::scalar Foam::kinematicParcelInjectionData::rho() const
{
    return rho_;
}


inline Foam::scalar Foam::kinematicParcelInjectionData::mDot() const
{
    return mDot_;
}


inline Foam::point& Foam::kinematicParcelInjectionData::x()
{
    return x_;
}


inline Foam::vector& Foam::kinematicParcelInjectionData::U()
{
    return U_;
}


inline Foam::scalar& Foam::kinematicParcelInjectionData::d()
{
    return d_;
}


inline Foam::scalar& Foam::kinematicParcelInjectionData::rho()
{
    return rho_;
}


inline Foam::scalar& Foam::kinematicParcelInjectionData::mDot()
{
    return mDot_;
}


// ************************************************************************* //
