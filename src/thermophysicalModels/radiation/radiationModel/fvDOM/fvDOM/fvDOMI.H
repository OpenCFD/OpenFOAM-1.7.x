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

inline const Foam::radiation::radiativeIntensityRay&
Foam::radiation::fvDOM::IRay(const label rayI) const
{
    return  IRay_[rayI];
}


inline const Foam::volScalarField&
Foam::radiation::fvDOM::IRayLambda
(
    const label rayI,
    const label lambdaI
) const
{
    return IRay_[rayI].ILambda(lambdaI);
}


inline Foam::label Foam::radiation::fvDOM::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::radiation::fvDOM::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::radiation::fvDOM::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::radiation::fvDOM::nLambda() const
{
    return nLambda_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOM::a() const
{
    return a_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOM::aLambda
(
    const label lambdaI
) const
{
    return aLambda_[lambdaI];
}


inline const Foam::volScalarField& Foam::radiation::fvDOM::G() const
{
    return G_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOM::Qr() const
{
    return Qr_;
}


inline const Foam::radiation::blackBodyEmission&
Foam::radiation::fvDOM::blackBody() const
{
    return blackBody_;
}


// ************************************************************************* //
