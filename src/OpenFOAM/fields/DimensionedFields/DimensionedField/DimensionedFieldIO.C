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

#include "DimensionedField.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type, class GeoMesh>
void DimensionedField<Type, GeoMesh>::readField
(
    const dictionary& fieldDict,
    const word& fieldDictEntry
)
{
    dimensions_.reset(dimensionSet(fieldDict.lookup("dimensions")));

    Field<Type> f(fieldDictEntry, fieldDict, GeoMesh::size(mesh_));
    transfer(f);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
DimensionedField<Type, GeoMesh>::DimensionedField
(
    const IOobject& io,
    const Mesh& mesh,
    const word& fieldDictEntry
)
:
    regIOobject(io),
    Field<Type>(0),
    mesh_(mesh),
    dimensions_(dimless)
{
    readField(dictionary(readStream(typeName)), fieldDictEntry);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMesh>
bool DimensionedField<Type, GeoMesh>::writeData
(
    Ostream& os,
    const word& fieldDictEntry
) const
{
    os.writeKeyword("dimensions") << dimensions() << token::END_STATEMENT
        << nl << nl;

    Field<Type>::writeEntry(fieldDictEntry, os);

    // Check state of Ostream
    os.check
    (
        "bool DimensionedField<Type, GeoMesh>::writeData"
        "(Ostream& os, const word& fieldDictEntry) const"
    );

    return (os.good());
}


template<class Type, class GeoMesh>
bool DimensionedField<Type, GeoMesh>::writeData(Ostream& os) const
{
    return writeData(os, "value");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type, class GeoMesh>
Ostream& operator<<(Ostream& os, const DimensionedField<Type, GeoMesh>& df)
{
    df.writeData(os);

    return os;
}


template<class Type, class GeoMesh>
Ostream& operator<<
(
    Ostream& os,
    const tmp<DimensionedField<Type, GeoMesh> >& tdf
)
{
    tdf().writeData(os);
    tdf.clear();

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
