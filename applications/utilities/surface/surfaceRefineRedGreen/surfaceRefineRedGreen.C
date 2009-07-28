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
    surfaceRefineRedGreen

Description
    Refine by splitting all three edges of triangle ('red' refinement).
    Neighbouring triangles (which are not marked for refinement get split
    in half ('green') refinement. (R. Verfuerth, "A review of a posteriori
    error estimation and adaptive mesh refinement techniques",
    Wiley-Teubner, 1996)

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "triSurfaceTools.H"
#include "argList.H"
#include "OFstream.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("surface file");
    argList::validArgs.append("output surface file");
    argList::argList args(argc, argv);

    fileName surfFileName(args.additionalArgs()[0]);
    fileName outFileName(args.additionalArgs()[1]);

    Info<< "Reading surface from " << surfFileName << " ..." << endl;

    triSurface surf1(surfFileName);

    // Refine
    triSurface surf2 = triSurfaceTools::redGreenRefine
    (
        surf1,
        identity(surf1.size())  //Hack: refine all
    );

    Info<< "Original surface:" << endl
        << "    triangles     :" << surf1.size() << endl
        << "    vertices(used):" << surf1.nPoints() << endl
        << "Refined surface:" << endl
        << "    triangles     :" << surf2.size() << endl
        << "    vertices(used):" << surf2.nPoints() << endl << endl;


    Info<< "Writing refined surface to " << outFileName << " ..." << endl;

    surf2.write(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
