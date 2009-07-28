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

\*---------------------------------------------------------------------------*/

#include "distributedTriSurfaceMesh.H"
#include "triSurfaceFields.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//template<class Type>
//void Foam::distributedTriSurfaceMesh::getField
//(
//    const word& fieldName,
//    const List<pointIndexHit>& info,
//    List<Type>& values
//) const
//{
//    typedef DimensionedField<Type, triSurfaceGeoMesh> DimensionedSurfField;
//
//
//    // Get query data (= local index of triangle)
//    // ~~~~~~~~~~~~~~
//
//    labelList triangleIndex(info.size());
//    autoPtr<mapDistribute> mapPtr
//    (
//        calcLocalQueries
//        (
//            info,
//            triangleIndex
//        )
//    );
//    const mapDistribute& map = mapPtr();
//
//
//    // Do my tests
//    // ~~~~~~~~~~~
//
//    const DimensionedSurfField& fld = lookupObject<DimensionedSurfField>
//    (
//        fieldName
//    );
//    const triSurface& s = static_cast<const triSurface&>(*this);
//
//    values.setSize(triangleIndex.size());
//
//    forAll(triangleIndex, i)
//    {
//        label triI = triangleIndex[i];
//        values[i] = fld[triI];
//    }
//
//
//    // Send back results
//    // ~~~~~~~~~~~~~~~~~
//
//    map.distribute
//    (
//        Pstream::nonBlocking,
//        List<labelPair>(0),
//        info.size(),
//        map.constructMap(),     // what to send
//        map.subMap(),           // what to receive
//        values
//    );
//}


template<class Type>
void Foam::distributedTriSurfaceMesh::distributeFields
(
    const mapDistribute& map
)
{
    typedef DimensionedField<Type, triSurfaceGeoMesh> DimensionedSurfField;

    HashTable<const DimensionedSurfField*> fields
    (
        objectRegistry::lookupClass
        <DimensionedSurfField >()
    );

    for
    (
        typename HashTable<const DimensionedSurfField*>::iterator fieldIter =
            fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        DimensionedSurfField& field =
            const_cast<DimensionedSurfField&>(*fieldIter());

        label oldSize = field.size();

        map.distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(0),
            map.constructSize(),
            map.subMap(),
            map.constructMap(),
            field
        );

        if (debug)
        {
            Info<< "Mapped " << field.typeName << ' ' << field.name()
                << " from size " << oldSize << " to size " << field.size()
                << endl;
        }
    }
}


// ************************************************************************* //
