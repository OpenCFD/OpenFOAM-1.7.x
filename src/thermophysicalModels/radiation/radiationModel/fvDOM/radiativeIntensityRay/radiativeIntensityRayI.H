/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

inline const Foam::volScalarField& Foam::radiation::radiativeIntensityRay::I() const
{
    return I_;
}


inline const Foam::volScalarField& Foam::radiation::radiativeIntensityRay::Qr() const
{
    return Qr_;
}


inline Foam::volScalarField& Foam::radiation::radiativeIntensityRay::Qr()
{
    return Qr_;
}


inline const Foam::vector& Foam::radiation::radiativeIntensityRay::d() const
{
    return d_;
}


inline const Foam::vector& Foam::radiation::radiativeIntensityRay::dAve() const
{
    return dAve_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRay::nLambda() const
{
    return nLambda_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRay::phi() const
{
    return phi_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRay::theta() const
{
    return theta_;
}


inline Foam::scalar Foam::radiation::radiativeIntensityRay::omega() const
{
    return omega_;
}


inline const Foam::volScalarField& Foam::radiation::radiativeIntensityRay::ILambda
(
    const label lambdaI
) const
{
    return ILambda_[lambdaI];
}


// ************************************************************************* //
