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

Application
    analyticalCylinder

Description
    Generates an analytical solution for potential flow around a cylinder.  
    Can be compared with the solution from the potentialFlow/cylinder example.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nEvaluating analytical solution" << endl;

    volVectorField centres = UA.mesh().C();
    volScalarField magCentres = mag(centres);
    volScalarField theta = acos((centres & vector(1,0,0))/magCentres);

    volVectorField cs2theta =
        cos(2*theta)*vector(1,0,0) 
      + sin(2*theta)*vector(0,1,0);

    UA = uInfX*(dimensionedVector(vector(1,0,0))
      - pow((radius/magCentres),2)*cs2theta);

    // Force writing of UA (since time has not changed)
    UA.write();

    Info<< "end" << endl;

    return 0;
}

// ************************************************************************* //
