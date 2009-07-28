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

Application
    graphTest

Description
    Test program for making graphs

\*---------------------------------------------------------------------------*/

#include "graph.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main()
{
    scalarField eta(200);
    scalarField C1mR(eta.size());

    forAll(eta, i)
    {
        eta[i] = scalar(i)/10.0;
    }

    scalar C1 = 1.42;
    scalar eta0 = 4.38;
    scalar beta = 0.012;

    C1mR = C1 - ((eta*(1.0 - eta/eta0))/(1.0 + beta*pow(eta, 3.0)));

    graph("C&1! - R", "C&1! - R", "@h!", eta,  C1mR).write
    (
        "C1mR",
        "xmgr"
    );

    Info<< "end" << endl;
}


// ************************************************************************* //

