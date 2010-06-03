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

Class
    Foam::ZoneMesh

Description
    List of mesh zones

SourceFiles
    ZoneMesh.C

\*---------------------------------------------------------------------------*/

#ifndef ZoneMesh_H
#define ZoneMesh_H

#include "List.H"
#include "regIOobject.H"
#include "HashSet.H"
#include "pointFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declaration of friend functions and operators

template<class ZoneType, class MeshType> class ZoneMesh;

template<class ZoneType, class MeshType>
Ostream& operator<<(Ostream&, const ZoneMesh<ZoneType, MeshType>&);

/*---------------------------------------------------------------------------*\
                           Class ZoneMesh Declaration
\*---------------------------------------------------------------------------*/

template<class ZoneType, class MeshType>
class ZoneMesh
:
    public PtrList<ZoneType>,
    public regIOobject
{
    // Private data

        //- Reference to mesh
        const MeshType& mesh_;

        //- Map of zone labels for given element
        mutable Map<label>* zoneMapPtr_;


    // Private member functions

        //- Disallow construct as copy
        ZoneMesh(const ZoneMesh&);

        //- Disallow assignment
        void operator=(const ZoneMesh<ZoneType, MeshType>&);


        //- Create zone map
        void calcZoneMap() const;


public:

    // Constructors

        //- Read constructor given IOobject and a MeshType reference
        ZoneMesh
        (
            const IOobject&,
            const MeshType&
        );

        //- Construct given size
        ZoneMesh
        (
            const IOobject&,
            const MeshType&,
            const label size
        );

    // Destructor

        ~ZoneMesh();


    // Member functions

        //- Return the mesh reference
        const MeshType& mesh() const
        {
            return mesh_;
        }

        //- Map of zones containing zone index for all zoned elements
        //  Return -1 if the object is not in the zone
        const Map<label>& zoneMap() const;

        //- Given a global object index, return the zone it is in.  If
        //object does not belong to any zones, return -1
        label whichZone(const label objectIndex) const;

        //- Return a list of zone types
        wordList types() const;

        //- Return a list of zone names
        wordList names() const;

        //- Find zone index given a name
        label findZoneID(const word& zoneName) const;

        //- Clear addressing
        void clearAddressing();

        //- Clear the zones
        void clear();

        //- Check zone definition. Return true if in error.
        bool checkDefinition(const bool report = false) const;

        //- Correct zone mesh after moving points
        void movePoints(const pointField&);

        //- writeData member function required by regIOobject
        bool writeData(Ostream&) const;


    // Ostream operator

        friend Ostream& operator<< <ZoneType, MeshType>
        (
            Ostream&,
            const ZoneMesh<ZoneType, MeshType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ZoneMesh.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
