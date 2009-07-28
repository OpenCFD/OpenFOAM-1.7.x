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
    foamFormatConvert

Description
    Converts all IOobjects associated with a case into the format specified
    in the controlDict.

    Mainly used to convert binary mesh/field files to ASCII.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "cellIOList.H"
#include "IOobjectList.H"
#include "IOPtrList.H"

#include "writeMeshObject.H"
#include "fieldDictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<entry>, 0);
}


// Main program:

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
#   include "setRootCase.H"
#   include "createTime.H"
    Foam::instantList timeDirs = Foam::timeSelector::select0(runTime, args);

    forAll (timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        // Convert all the standard mesh files
        writeMeshObject<cellIOList>("cells", runTime);
        writeMeshObject<labelIOList>("owner", runTime);
        writeMeshObject<labelIOList>("neighbour", runTime);
        writeMeshObject<faceIOList>("faces", runTime);
        writeMeshObject<pointIOField>("points", runTime);
        writeMeshObject<IOPtrList<entry> >("cellZones", runTime);
        writeMeshObject<IOPtrList<entry> >("faceZones", runTime);
        writeMeshObject<IOPtrList<entry> >("pointZones", runTime);

        // Get list of objects from the database
        IOobjectList objects(runTime, runTime.timeName());

        forAllConstIter(IOobjectList, objects, iter)
        {
            const word& headerClassName = iter()->headerClassName();

            if
            (
                headerClassName == volScalarField::typeName
             || headerClassName == volVectorField::typeName
             || headerClassName == volSphericalTensorField::typeName
             || headerClassName == volSymmTensorField::typeName
             || headerClassName == volTensorField::typeName

             || headerClassName == surfaceScalarField::typeName
             || headerClassName == surfaceVectorField::typeName
             || headerClassName == surfaceSphericalTensorField::typeName
             || headerClassName == surfaceSymmTensorField::typeName
             || headerClassName == surfaceTensorField::typeName

             || headerClassName == pointScalarField::typeName
             || headerClassName == pointVectorField::typeName
             || headerClassName == pointSphericalTensorField::typeName
             || headerClassName == pointSymmTensorField::typeName
             || headerClassName == pointTensorField::typeName
            )
            {
                Info<< "        Reading " << headerClassName
                    << " : " << iter()->name() << endl;

                fieldDictionary fDict
                (
                    *iter(),
                    headerClassName
                );

                Info<< "        Writing " << iter()->name() << endl;
                fDict.regIOobject::write();
            }
        }

        Info<< endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
