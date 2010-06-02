/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Description
    Translates FOAM data to Fluent format.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "writeFluentFields.H"
#include "OFstream.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions(false);   // no constant

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    // make a directory called proInterface in the case
    mkDir(runTime.rootPath()/runTime.caseName()/"fluentInterface");

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        if (mesh.readUpdate())
        {
            Info<< "    Read new mesh" << endl;
        }

        // make a directory called proInterface in the case
        mkDir(runTime.rootPath()/runTime.caseName()/"fluentInterface");

        // open a file for the mesh
        OFstream fluentDataFile
        (
            runTime.rootPath()/
            runTime.caseName()/
            "fluentInterface"/
            runTime.caseName() + runTime.timeName() + ".dat"
        );

        fluentDataFile
            << "(0 \"FOAM to Fluent data File\")" << endl << endl;

        // Writing number of faces
        label nFaces = mesh.nFaces();

        forAll (mesh.boundary(), patchI)
        {
            nFaces += mesh.boundary()[patchI].size();
        }

        fluentDataFile
            << "(33 (" << mesh.nCells() << " " << nFaces << " "
            << mesh.nPoints() << "))" << endl;

        IOdictionary foamDataToFluentDict
        (
            IOobject
            (
                "foamDataToFluentDict",
                runTime.system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );


        // Search for list of objects for this time
        IOobjectList objects(mesh, runTime.timeName());


        // Converting volScalarField
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Search list of objects for volScalarFields
        IOobjectList scalarFields(objects.lookupClass("volScalarField"));

        for
        (
            IOobjectList::iterator scalarFieldIter = scalarFields.begin();
            scalarFieldIter != scalarFields.end();
            ++scalarFieldIter
        )
        {
            // Read field
            volScalarField field
            (
                *scalarFieldIter(),
                mesh
            );

            // lookup field from dictionary
            if (foamDataToFluentDict.found(field.name()))
            {
                label unitNumber
                (
                    readLabel(foamDataToFluentDict.lookup(field.name()))
                );

                // Convert field
                if (unitNumber > 0)
                {
                    Info<< "    Converting field " << field.name() << endl;
                    writeFluentField(field, unitNumber, fluentDataFile);
                }
            }
        }


        // Converting volVectorField
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Search list of objects for volVectorFields
        IOobjectList vectorFields(objects.lookupClass("volVectorField"));

        for
        (
            IOobjectList::iterator vectorFieldIter = vectorFields.begin();
            vectorFieldIter != vectorFields.end();
            ++vectorFieldIter
        )
        {
            // Read field
            volVectorField field
            (
                *vectorFieldIter(),
                mesh
            );

            // lookup field from dictionary
            if (foamDataToFluentDict.found(field.name()))
            {
                label unitNumber
                (
                    readLabel(foamDataToFluentDict.lookup(field.name()))
                );

                // Convert field
                if (unitNumber > 0)
                {
                    Info<< "    Converting field " << field.name() << endl;
                    writeFluentField(field, unitNumber, fluentDataFile);
                }
            }
        }

        Info<< endl;
    }

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
