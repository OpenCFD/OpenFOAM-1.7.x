/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    testDirectMappedPatch

Description
    Test direct mapped b.c. by mapping face centres (mesh.C().boundaryField()).

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "meshTools.H"
#include "Time.H"
#include "OFstream.H"
#include "volFields.H"
#include "directMappedFixedValueFvPatchFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Main program:

int main(int argc, char *argv[])
{
#   include "addTimeOptions.H"
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    wordList patchFieldTypes
    (
        mesh.boundaryMesh().size(),
        calculatedFvPatchVectorField::typeName
    );

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<directMappedPolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            patchFieldTypes[patchI] =
                directMappedFixedValueFvPatchVectorField::typeName;
        }
    }

    Pout<< "patchFieldTypes:" << patchFieldTypes << endl;

    volVectorField cc
    (
        IOobject
        (
            "cc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero),
        patchFieldTypes
    );
    
    cc.internalField() = mesh.C().internalField();
    cc.boundaryField().updateCoeffs();

    forAll(cc.boundaryField(), patchI)
    {
        if
        (
            isA<directMappedFixedValueFvPatchVectorField>
            (
                cc.boundaryField()[patchI]
            )
        )
        {
            Pout<< "Detected a directMapped patch:" << patchI << endl;

            OFstream str(mesh.boundaryMesh()[patchI].name() + ".obj");
            Pout<< "Writing mapped values to " << str.name() << endl;

            label vertI = 0;
            const fvPatchVectorField& fvp = cc.boundaryField()[patchI];

            forAll(fvp, i)
            {
                meshTools::writeOBJ(str, fvp.patch().Cf()[i]);
                vertI++;
                meshTools::writeOBJ(str, fvp[i]);
                vertI++;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
