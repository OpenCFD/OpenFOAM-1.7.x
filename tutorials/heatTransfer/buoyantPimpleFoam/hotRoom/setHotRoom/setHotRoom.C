/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    setHotRoom

Description
    Set the initial field of T for the hot room problem.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OSspecific.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

volScalarField::GeometricBoundaryField& Tpatches = T.boundaryField();

forAll(Tpatches, patchI)
{
    if
    (
        isA<fixedValueFvPatchScalarField>(Tpatches[patchI])
     && mesh.boundaryMesh()[patchI].name() == "floor"
    )
    {
        fixedValueFvPatchScalarField& Tpatch =
            refCast<fixedValueFvPatchScalarField>(Tpatches[patchI]);
        
        const vectorField& faceCentres = 
            mesh.Cf().boundaryField()[patchI];

        forAll(faceCentres, facei)
        {
            if 
            (
                (faceCentres[facei].x() > 4.5) && 
                (faceCentres[facei].x() < 5.5) && 
                (faceCentres[facei].z() > 4.5) && 
                (faceCentres[facei].z() < 5.5) 
            )
            {
                Tpatch[facei] = 600;
            }
            else
            {
                Tpatch[facei] = 300;
            }
        }
    };

    Info<< "Writing modified field T\n" << endl;
    T.write();

    Info<< "End\n" << endl;

    return 0;
}

}
// ************************************************************************* //
