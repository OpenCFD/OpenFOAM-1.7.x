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

Application
    changeDictionary

Description
    Utility to change dictionary entries, e.g. can be used to change the patch
    type in the field and polyMesh/boundary files.

    Reads dictionaries (fields) and entries to change from a dictionary.
    E.g. to make the @em movingWall a @em fixedValue for @em p but all other
    @em Walls a zeroGradient boundary condition, the
    @c system/changeDictionaryDict would contain the following:
    @verbatim
    dictionaryReplacement
    {
        p                           // field to change
        {
            boundaryField
            {
                ".*Wall"            // entry to change
                {
                    type            zeroGradient;
                }
                movingWall          // entry to change
                {
                    type            fixedValue;
                    value           uniform 123.45;
                }
            }
        }
    }
    @endverbatim

Usage

    - changeDictionary [OPTION]

    @param -literalRE \n
    Do not interpret regular expressions; treat them as any other keyword.


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IOobjectList.H"
#include "IOPtrList.H"
#include "volFields.H"
#include "stringListOps.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(IOPtrList<entry>, 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool merge(dictionary&, const dictionary&, const bool);


// Add thisEntry to dictionary thisDict.
bool addEntry
(
    dictionary& thisDict,
    entry& thisEntry,
    const entry& mergeEntry,
    const bool literalRE
)
{
    bool changed = false;

    // Recursively merge sub-dictionaries
    // TODO: merge without copying
    if (thisEntry.isDict() && mergeEntry.isDict())
    {
        if
        (
            merge
            (
                const_cast<dictionary&>(thisEntry.dict()),
                mergeEntry.dict(),
                literalRE
            )
        )
        {
            changed = true;
        }
    }
    else
    {
        // Should use in-place modification instead of adding
        thisDict.add(mergeEntry.clone(thisDict).ptr(), true);
        changed = true;
    }

    return changed;
}


// Dictionary merging/editing.
// literalRE:
// - true: behave like dictionary::merge, i.e. add regexps just like
//   any other key.
// - false : interpret wildcard as a rule for items to be matched.
bool merge
(
    dictionary& thisDict,
    const dictionary& mergeDict,
    const bool literalRE
)
{
    auto bool wildCardInMergeDict = false;

    bool changed = false;

    // Save current (non-wildcard) keys before adding items.
    HashSet<word> thisKeysSet;
    {
        List<keyType> keys = thisDict.keys(false);
        forAll(keys, i)
        {
            thisKeysSet.insert(keys[i]);
        }
    }

    // Pass 1. All literal matches

    forAllConstIter(IDLList<entry>, mergeDict, mergeIter)
    {
        const keyType& key = mergeIter().keyword();

        if (literalRE || !key.isPattern())
        {
            entry* entryPtr = thisDict.lookupEntryPtr
            (
                key,
                false,              // recursive
                false               // patternMatch
            );

            if (entryPtr)
            {

                // Mark thisDict entry as having been match for wildcard
                // handling later on.
                thisKeysSet.erase(entryPtr->keyword());

                if
                (
                    addEntry
                    (
                        thisDict,
                       *entryPtr,
                        mergeIter(),
                        literalRE
                    )
                )
                {
                    changed = true;
                }
            }
            else
            {
                // not found - just add
                thisDict.add(mergeIter().clone(thisDict).ptr());
                changed = true;
            }
        }
    }


    // Pass 2. Wildcard matches (if any) on any non-match keys.

    if (!literalRE && thisKeysSet.size() > 0)
    {
        wordList thisKeys(thisKeysSet.toc());

        forAllConstIter(IDLList<entry>, mergeDict, mergeIter)
        {
            const keyType& key = mergeIter().keyword();

            if (key.isPattern())
            {
                // Find all matching entries in the original thisDict

                if (!wildCardInMergeDict)
                {
                    wildCardInMergeDict = true;
                    WarningIn("changeDictionary()")
                        << "Detected wildcard " << key
                        << " in changeDictionaryDict" << endl
                        << "The behaviour of wildcards has changed -"
                        << " they are now interpreted by changeDictionary."
                        << endl << "Please take care or use the -literalRE"
                        << " command line option to revert to"
                        << " previous behaviour." << endl;
                }

                labelList matches = findStrings(key, thisKeys);

                forAll(matches, i)
                {
                    label matchI = matches[i];

                    entry& thisEntry = const_cast<entry&>
                    (
                        thisDict.lookupEntry(thisKeys[matchI], false, false)
                    );

                    if
                    (
                        addEntry
                        (
                            thisDict,
                            thisEntry,
                            mergeIter(),
                            literalRE
                        )
                    )
                    {
                        changed = true;
                    }
                }
            }
        }
    }

    return changed;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("instance", "instance");
    argList::validOptions.insert("literalRE", "");
    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    bool literalRE = args.optionFound("literalRE");

    if (literalRE)
    {
        Info<< "Not interpreting any regular expressions (RE)"
            << " in the changeDictionaryDict." << endl
            << "Instead they are handled as any other entry, i.e. added if"
            << " not present." << endl;
    }


    fileName regionPrefix = "";
    if (regionName != fvMesh::defaultRegion)
    {
        regionPrefix = regionName;
    }

    word instance = runTime.timeName();
    if (args.options().found("instance"))
    {
        instance = args.options()["instance"];
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
            merge(fieldDict, replaceDict, literalRE);

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
                    instance,
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
            merge(fieldDict, replaceDict, literalRE);

            Info<< "Writing modified fieldDict " << fieldName << endl;
            fieldDict.regIOobject::write();
        }
    }

    Info<< endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
