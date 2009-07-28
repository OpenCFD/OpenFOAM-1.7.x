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

\*---------------------------------------------------------------------------*/

#include "compressibleCourantNo.H"
#include "fvc.H"

Foam::scalar Foam::compressibleCourantNo
(
    const fvMesh& mesh,
    const Time& runTime,
    const volScalarField& rho,
    const surfaceScalarField& phi
)
{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    //- Can have fluid domains with 0 cells so do not test.
    //if (mesh.nInternalFaces())
    {
        surfaceScalarField SfUfbyDelta =
            mesh.surfaceInterpolation::deltaCoeffs()
          * mag(phi)
          / fvc::interpolate(rho);

        CoNum = max(SfUfbyDelta/mesh.magSf())
            .value()*runTime.deltaT().value();

        meanCoNum = (sum(SfUfbyDelta)/sum(mesh.magSf()))
            .value()*runTime.deltaT().value();
    }

    Info<< "Region: " << mesh.name() << " Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;

    return CoNum;
}


// ************************************************************************* //
