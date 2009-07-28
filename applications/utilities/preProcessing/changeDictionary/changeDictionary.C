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
    changeDictionary

Description
    Utility to change dictionary entries, e.g. can be used to change the patch
    type in the field and polyMesh/boundary files.

    Reads dictionaries (fields) and entries to change from a dictionary.
    E.g. to make the @em movingWall a @em fixedValue for @em p, the
    @c system/changeDictionaryDict would contain the following:
    @verbatim
    dictionaryReplacement
    {
        p                           // field to change
        {
            boundaryField
            {
                movingWall          // entry to change
                {
                    type            fixedValue;
                    value           uniform 123.45;
                }
            }
        }
    }
    @endverbatim


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "volFields.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<entry>, 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "addRegionOption.H"
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createNamedMesh.H"

    fileName regionPrefix = "";
    if (regionName != fvMesh::defaultRegion)
    {
        regionPrefix = regionName;
    }

    // Get the replacement rules from a dictionary
    IOdictionary dict
    (
        IOobject
        (
            "changeDictionaryDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    const dictionary& replaceDicts = dict.subDict("dictionaryReplacement");
    Info<< "Read dictionary " << dict.name()
        << " with replacements for dictionaries "
        << replaceDicts.toc() << endl;


    // Every replacement is a dictionary name and a keyword in this

    forAllConstIter(dictionary, replaceDicts, fieldIter)
    {
        const word& fieldName = fieldIter().keyword();
        Info<< "Replacing entries in dictionary " << fieldName << endl;

        // Handle 'boundary' specially:
        // - is PtrList of dictionaries
        // - is in polyMesh/
        if (fieldName == "boundary")
        {
            Info<< "Special handling of " << fieldName
                << " as polyMesh/boundary file." << endl;

            // Read PtrList of dictionary as dictionary.
            const word oldTypeName = IOPtrList<entry>::typeName;
            const_cast<word&>(IOPtrList<entry>::typeName) = word::null;
            IOPtrList<entry> dictList
            (
                IOobject
                (
                    fieldName,
                    runTime.findInstance
                    (
                        regionPrefix/polyMesh::meshSubDir,
                        fieldName
                    ),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );
            const_cast<word&>(IOPtrList<entry>::typeName) = oldTypeName;
            // Fake type back to what was in field
            const_cast<word&>(dictList.type()) = dictList.headerClassName();

            // Temporary convert to dictionary
            dictionary fieldDict;
            forAll(dictList, i)
            {
                fieldDict.add(dictList[i].keyword(), dictList[i].dict());
            }

            Info<< "Loaded dictionary " << fieldName
                << " with entries " << fieldDict.toc() << endl;

            // Get the replacement dictionary for the field
            const dictionary& replaceDict = fieldIter().dict();
            Info<< "Merging entries from " << replaceDict.toc() << endl;

            // Merge the replacements in
            fieldDict.merge(replaceDict);

            Info<< "fieldDict:" << fieldDict << endl;

            // Convert back into dictList
            wordList doneKeys(dictList.size());

            label nEntries = fieldDict.size();
            forAll(dictList, i)
            {
                doneKeys[i] = dictList[i].keyword();
                dictList.set
                (
                    i,
                    fieldDict.lookupEntry
                    (
                        doneKeys[i],
                        false,
                        true
                    ).clone()
                );
                fieldDict.remove(doneKeys[i]);
            }
            // Add remaining entries
            label sz = dictList.size();
            dictList.setSize(nEntries);
            forAllConstIter(dictionary, fieldDict, iter)
            {
                dictList.set(sz, iter().clone());
            }

            Info<< "Writing modified fieldDict " << fieldName << endl;
            dictList.write();
        }
        else
        {
            // Read dictionary. (disable class type checking so we can load
            // field)
            Info<< "Loading dictionary " << fieldName << endl;
            const word oldTypeName = IOdictionary::typeName;
            const_cast<word&>(IOdictionary::typeName) = word::null;
            IOdictionary fieldDict
            (
                IOobject
                (
                    fieldName,
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );
            const_cast<word&>(IOdictionary::typeName) = oldTypeName;
            // Fake type back to what was in field
            const_cast<word&>(fieldDict.type()) = fieldDict.headerClassName();

            Info<< "Loaded dictionary " << fieldName
                << " with entries " << fieldDict.toc() << endl;

            // Get the replacement dictionary for the field
            const dictionary& replaceDict = fieldIter().dict();
            Info<< "Merging entries from " << replaceDict.toc() << endl;

            // Merge the replacements in
            fieldDict.merge(replaceDict);

            Info<< "Writing modified fieldDict " << fieldName << endl;
            fieldDict.regIOobject::write();
        }
    }

    Info<< endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
