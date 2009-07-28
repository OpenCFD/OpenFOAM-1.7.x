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
    patchSummary

Description
    Writes fields and boundary condition info for each patch at each requested
    time instance.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "volFields.H"
#include "IOobjectList.H"
#include "patchSummaryTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << nl << endl;

        const IOobjectList fieldObjs(mesh, runTime.timeName());
        const wordList objNames = fieldObjs.names();

        PtrList<volScalarField> vsf(objNames.size());
        PtrList<volVectorField> vvf(objNames.size());
        PtrList<volSphericalTensorField> vsptf(objNames.size());
        PtrList<volSymmTensorField> vsytf(objNames.size());
        PtrList<volTensorField> vtf(objNames.size());

        Info<< "Valid fields:" << endl;

        forAll(objNames, objI)
        {
            IOobject obj
            (
                objNames[objI],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ
            );

            if (obj.headerOk())
            {
                addToFieldList<scalar>(vsf, obj, objI, mesh);
                addToFieldList<vector>(vvf, obj, objI, mesh);
                addToFieldList<sphericalTensor>(vsptf, obj, objI, mesh);
                addToFieldList<symmTensor>(vsytf, obj, objI, mesh);
                addToFieldList<tensor>(vtf, obj, objI, mesh);
            }
        }

        Info<< endl;

        const polyBoundaryMesh& bm = mesh.boundaryMesh();
        forAll(bm, patchI)
        {
            Info<< bm[patchI].type() << ": " << bm[patchI].name() << nl;
            outputFieldList<scalar>(vsf, patchI);
            outputFieldList<vector>(vvf, patchI);
            outputFieldList<sphericalTensor>(vsptf, patchI);
            outputFieldList<symmTensor>(vsytf, patchI);
            outputFieldList<tensor>(vtf, patchI);
            Info << endl;
        }
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
