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
    Translates FOAM data to EnSight format.

    An Ensight part is created for the internalMesh and for each patch.

Usage
    - foamToEnsight [OPTION] \n
    Translates OpenFOAM data to Ensight format

    @param -ascii \n
    Write Ensight data in ASCII format instead of "C Binary"

    @param -patches patchList \n
    Specify particular patches to write.
    Specifying an empty list suppresses writing the internalMesh.

    @param -noPatches \n
    Suppress writing any patches.

Note
    Parallel support for cloud data is not supported
    - writes to @a EnSight directory to avoid collisions with foamToEnsightParts

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "IOmanip.H"
#include "OFstream.H"

#include "volFields.H"

#include "labelIOField.H"
#include "scalarIOField.H"
#include "tensorIOField.H"

#include "ensightMesh.H"
#include "ensightField.H"

#include "ensightParticlePositions.H"
#include "ensightCloudField.H"

#include "fvc.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool inFileNameList
(
    const fileNameList& nameList,
    const word& name
)
{
    forAll(nameList, i)
    {
        if (nameList[i] == name)
        {
            return true;
        }
    }

    return false;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("ascii", "" );
    argList::validOptions.insert("patches", "patchList");
    argList::validOptions.insert("noPatches", "");

#   include "addTimeOptions.H"
#   include "setRootCase.H"

    // Check options
    bool binary = !args.optionFound("ascii");

#   include "createTime.H"

    // get the available time-steps
    instantList Times = runTime.times();

#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createNamedMesh.H"

    // Mesh instance (region0 gets filtered out)
    fileName regionPrefix = "";

    if (regionName != polyMesh::defaultRegion)
    {
        regionPrefix = regionName;
    }

    const label nVolFieldTypes = 5;
    const word volFieldTypes[] =
    {
        volScalarField::typeName,
        volVectorField::typeName,
        volSphericalTensorField::typeName,
        volSymmTensorField::typeName,
        volTensorField::typeName
    };

    // Path to EnSight folder at case level only
    // - For parallel cases, data only written from master
    fileName ensightDir = args.rootPath()/args.globalCaseName()/"EnSight";

    if (Pstream::master())
    {
        if (isDir(ensightDir))
        {
            rmDir(ensightDir);
        }

        mkDir(ensightDir);
    }

    // Start of case file header output
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const word prepend = args.globalCaseName() + '.';

    OFstream *ensightCaseFilePtr = NULL;
    if (Pstream::master())
    {
        fileName caseFileName = prepend + "case";
        Info<< nl << "write case: " << caseFileName.c_str() << endl;

        // the case file is always ASCII
        ensightCaseFilePtr = new OFstream
        (
            ensightDir/caseFileName,
            IOstream::ASCII
        );

        *ensightCaseFilePtr
            << "FORMAT" << nl
            << "type: ensight gold" << nl << nl;
    }

    OFstream& ensightCaseFile = *ensightCaseFilePtr;

    // Construct the EnSight mesh
    ensightMesh eMesh(mesh, args, binary);

    // Set Time to the last time before looking for the lagrangian objects
    runTime.setTime(Times[Times.size()-1], Times.size()-1);

    IOobjectList objects(mesh, runTime.timeName());

#   include "checkMeshMoving.H"

    wordHashSet allCloudNames;
    if (Pstream::master())
    {
        word geomFileName = prepend + "000";

        // test pre check variable if there is a moving mesh
        if (meshMoving)
        {
            geomFileName = prepend + "***";
        }

        ensightCaseFile
            << "GEOMETRY" << nl
            << "model:        1     "
            << (geomFileName + ".mesh").c_str() << nl;
    }

    // Identify if lagrangian data exists at each time, and add clouds
    // to the 'allCloudNames' hash set
    for (label n=startTime; n<endTime; n++)
    {
        runTime.setTime(Times[n], n);

        fileNameList cloudDirs = readDir
        (
            runTime.timePath()/regionPrefix/cloud::prefix,
            fileName::DIRECTORY
        );

        forAll(cloudDirs, cloudI)
        {
            IOobjectList cloudObjs
            (
                mesh,
                runTime.timeName(),
                cloud::prefix/cloudDirs[cloudI]
            );

            IOobject* positionsPtr = cloudObjs.lookup("positions");

            if (positionsPtr)
            {
                allCloudNames.insert(cloudDirs[cloudI]);
            }
        }
    }

    HashTable<HashTable<word> > allCloudFields;
    forAllConstIter(wordHashSet, allCloudNames, cloudIter)
    {
        // Add the name of the cloud(s) to the case file header
        if (Pstream::master())
        {
            ensightCaseFile
            <<  (
                    "measured:     1     "
                  + prepend
                  + "***."
                  + cloudIter.key()
                ).c_str()
            << nl;
        }

        // Create a new hash table for each cloud
        allCloudFields.insert(cloudIter.key(), HashTable<word>());

        // Identify the new cloud in the hash table
        HashTable<HashTable<word> >::iterator newCloudIter =
            allCloudFields.find(cloudIter.key());

        // Loop over all times to build list of fields and field types
        // for each cloud
        for (label n=startTime; n<endTime; n++)
        {
            runTime.setTime(Times[n], n);

            IOobjectList cloudObjs
            (
                mesh,
                runTime.timeName(),
                cloud::prefix/cloudIter.key()
            );

            forAllConstIter(IOobjectList, cloudObjs, fieldIter)
            {
                const IOobject obj = *fieldIter();

                if (obj.name() != "positions")
                {
                    // Add field and field type
                    newCloudIter().insert
                    (
                        obj.name(),
                        obj.headerClassName()
                    );
                }
            }
        }
    }

    label nTimeSteps = 0;
    for (label n=startTime; n<endTime; n++)
    {
        nTimeSteps++;
        runTime.setTime(Times[n], n);
        label timeIndex = n - startTime;

        word timeName = itoa(timeIndex);
        word timeFile = prepend + timeName;

        Info<< "Translating time = " << runTime.timeName() << nl;

#       include "moveMesh.H"

        if (timeIndex == 0 || mesh.moving())
        {
            eMesh.write
            (
                ensightDir,
                prepend,
                timeIndex,
                ensightCaseFile
            );
        }


        // Start of field data output
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (timeIndex == 0 && Pstream::master())
        {
            ensightCaseFile<< nl << "VARIABLE" << nl;
        }


        // Cell field data output
        // ~~~~~~~~~~~~~~~~~~~~~~

        for (label i=0; i<nVolFieldTypes; i++)
        {
            wordList fieldNames = objects.names(volFieldTypes[i]);

            for (label j=0; j<fieldNames.size(); j++)
            {
                word fieldName = fieldNames[j];

#               include "checkData.H"

                if (!variableGood)
                {
                    continue;
                }

                IOobject fieldObject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                );

                if (volFieldTypes[i] == volScalarField::typeName)
                {
                    ensightField<scalar>
                    (
                        fieldObject,
                        eMesh,
                        ensightDir,
                        prepend,
                        timeIndex,
                        binary,
                        ensightCaseFile
                    );
                }
                else if (volFieldTypes[i] == volVectorField::typeName)
                {
                    ensightField<vector>
                    (
                        fieldObject,
                        eMesh,
                        ensightDir,
                        prepend,
                        timeIndex,
                        binary,
                        ensightCaseFile
                    );
                }
                else if (volFieldTypes[i] == volSphericalTensorField::typeName)
                {
                    ensightField<sphericalTensor>
                    (
                        fieldObject,
                        eMesh,
                        ensightDir,
                        prepend,
                        timeIndex,
                        binary,
                        ensightCaseFile
                    );
                }
                else if (volFieldTypes[i] == volSymmTensorField::typeName)
                {
                    ensightField<symmTensor>
                    (
                        fieldObject,
                        eMesh,
                        ensightDir,
                        prepend,
                        timeIndex,
                        binary,
                        ensightCaseFile
                    );
                }
                else if (volFieldTypes[i] == volTensorField::typeName)
                {
                    ensightField<tensor>
                    (
                        fieldObject,
                        eMesh,
                        ensightDir,
                        prepend,
                        timeIndex,
                        binary,
                        ensightCaseFile
                    );
                }
            }
        }


        // Cloud field data output
        // ~~~~~~~~~~~~~~~~~~~~~~~

        forAllConstIter(HashTable<HashTable<word> >, allCloudFields, cloudIter)
        {
            const word& cloudName = cloudIter.key();

            fileNameList currentCloudDirs = readDir
            (
                runTime.timePath()/regionPrefix/cloud::prefix,
                fileName::DIRECTORY
            );

            bool cloudExists = inFileNameList(currentCloudDirs, cloudName);
            ensightParticlePositions
            (
                mesh,
                ensightDir,
                timeFile,
                cloudName,
                cloudExists
            );

            forAllConstIter(HashTable<word>, cloudIter(), fieldIter)
            {
                const word& fieldName = fieldIter.key();
                const word& fieldType = fieldIter();

                IOobject fieldObject
                (
                    fieldName,
                    mesh.time().timeName(),
                    cloud::prefix/cloudName,
                    mesh,
                    IOobject::MUST_READ
                );

                bool fieldExists = fieldObject.headerOk();
                if (fieldType == scalarIOField::typeName)
                {
                    ensightCloudField<scalar>
                    (
                        fieldObject,
                        ensightDir,
                        prepend,
                        timeIndex,
                        cloudName,
                        ensightCaseFile,
                        fieldExists
                    );
                }
                else if (fieldType == vectorIOField::typeName)
                {
                    ensightCloudField<vector>
                    (
                        fieldObject,
                        ensightDir,
                        prepend,
                        timeIndex,
                        cloudName,
                        ensightCaseFile,
                        fieldExists
                    );
                }
                else
                {
                    Info<< "Unable to convert field type " << fieldType
                        << " for field " << fieldName << endl;
                }
            }
        }
    }

#   include "ensightCaseTail.H"

    if (Pstream::master())
    {
        delete ensightCaseFilePtr;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
