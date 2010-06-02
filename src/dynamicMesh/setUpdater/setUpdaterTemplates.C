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

\*---------------------------------------------------------------------------*/

#include "setUpdater.H"
#include "polyMesh.H"
#include "Time.H"
#include "mapPolyMesh.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void setUpdater::updateSets(const mapPolyMesh& morphMap) const
{
    //
    // Update all sets in memory.
    //

    HashTable<const Type*> memSets = 
        morphMap.mesh().objectRegistry::lookupClass<Type>();

    for
    (
        typename HashTable<const Type*>::iterator iter = memSets.begin();
        iter != memSets.end();
        ++iter
    )
    {
        Type& set = const_cast<Type&>(*iter());

        if (debug)
        {
            Pout<< "Set:" << set.name() << " size:" << set.size()
                << " updated in memory" << endl;
        }

        set.updateMesh(morphMap);

        // Write or not? Debatable.
        set.write();
    }


    //
    // Update all sets on disk
    //

    // Get last valid mesh (discard points-only change)
    IOobjectList Objects
    (
        morphMap.mesh().time(),
        morphMap.mesh().time().findInstance
        (
            morphMap.mesh().meshDir(),
            "faces"
        ),
        "polyMesh/sets"
    );

    IOobjectList fileSets(Objects.lookupClass(Type::typeName));

    for
    (
        IOobjectList::const_iterator iter = fileSets.begin();
        iter != fileSets.end();
        ++iter
    )
    {
        if (!memSets.found(iter.key()))
        {
            // Not in memory. Load it.
            Type set(*iter());

            if (debug)
            {
                Pout<< "Set:" << set.name() << " size:" << set.size()
                    << " updated on disk" << endl;
            }

            set.updateMesh(morphMap);

            set.write();
        }
        else
        {
            if (debug)
            {
                Pout<< "Set:" << iter.key() << " already updated from memory"
                    << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
