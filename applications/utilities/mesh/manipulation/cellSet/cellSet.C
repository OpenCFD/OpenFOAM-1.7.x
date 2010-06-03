/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Description
    Selects a cell set through a dictionary.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "topoSetSource.H"
#include "cellSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Copy set
void backup
(
    const polyMesh& mesh,
    const word& fromName,
    const topoSet& fromSet,
    const word& toName
)
{
    Info<< "Backing up " << fromName << " into " << toName << endl;

    topoSet backupSet(mesh, toName, fromSet);

    backupSet.write();
}


// Read and copy set
void backup
(
    const polyMesh& mesh,
    const word& fromName,
    const word& toName
)
{
    topoSet fromSet(mesh, fromName, IOobject::READ_IF_PRESENT);

    backup(mesh, fromName, fromSet, toName);
}


// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createPolyMesh.H"

    Info<< "Reading cellSetDict\n" << endl;

    IOdictionary cellSetDict
    (
        IOobject
        (
            "cellSetDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    word setName(cellSetDict.lookup("name"));

    word actionName(cellSetDict.lookup("action"));

    topoSetSource::setAction action = topoSetSource::toAction(actionName);


    // Create topoSetSources
    PtrList<topoSetSource> topoSetSources
    (
        cellSetDict.lookup("topoSetSources"),
        topoSetSource::iNew(mesh)
    );


    // Load set to work
    autoPtr<topoSet> currentSetPtr(NULL);
    IOobject::readOption r;

    if ((action == topoSetSource::NEW) || (action == topoSetSource::CLEAR))
    {
        r = IOobject::NO_READ;

        backup(mesh, setName, setName + "_old");

        currentSetPtr.reset
        (
            new cellSet
            (
                mesh,
                setName,
                mesh.nCells()/10+1  // Reasonable size estimate.
            )
        );
    }
    else
    {
        r = IOobject::MUST_READ;

        currentSetPtr.reset
        (
            new cellSet
            (
                mesh,
                setName,
                r
            )
        );
    }

    topoSet& currentSet = currentSetPtr();

    Info<< "Set:" << currentSet.name()
        << "  Size:" << currentSet.size()
        << "  Action:" << actionName
        << endl;

    if ((r == IOobject::MUST_READ) && (action != topoSetSource::LIST))
    {
        // currentSet has been read so can make copy.
        backup(mesh, setName, currentSet, setName + "_old");
    }

    if (action == topoSetSource::CLEAR)
    {
        // Already handled above by not reading
    }
    else if (action == topoSetSource::INVERT)
    {
        currentSet.invert(currentSet.maxSize(mesh));
    }
    else if (action == topoSetSource::LIST)
    {
        currentSet.writeDebug(Info, mesh, 100);
        Info<< endl;
    }
    else if (action == topoSetSource::SUBSET)
    {
        // Apply topoSetSources to it to handle new/add/delete
        forAll(topoSetSources, topoSetSourceI)
        {
            // Backup current set.
            topoSet oldSet(mesh, currentSet.name() + "_old2", currentSet);

            currentSet.clear();

            topoSetSources[topoSetSourceI].applyToSet
            (
                topoSetSource::NEW,
                currentSet
            );

            // Combine new value of currentSet with old one.
            currentSet.subset(oldSet);
        }
    }
    else
    {
        // Apply topoSetSources to it to handle new/add/delete
        forAll(topoSetSources, topoSetSourceI)
        {
            topoSetSources[topoSetSourceI].applyToSet(action, currentSet);
        }
    }


    if (action != topoSetSource::LIST)
    {
        // Set has changed.

        // Sync across coupled patches.
        currentSet.sync(mesh);

        Info<< "Writing " << currentSet.name()
            << " (size " << currentSet.size() << ") to "
            << currentSet.instance()/currentSet.local()
               /currentSet.name()
            << endl << endl;

        currentSet.write();
    }

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
