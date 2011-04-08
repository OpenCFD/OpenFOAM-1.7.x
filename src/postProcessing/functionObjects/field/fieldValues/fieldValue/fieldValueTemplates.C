/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "fieldValue.H"
#include "ListListOps.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fieldValue::combineFields
(
    const Field<Type>& field
) const
{
    List<Field<Type> > allValues(Pstream::nProcs());

    allValues[Pstream::myProcNo()] = field;

    Pstream::gatherList(allValues);

    if (Pstream::master())
    {
        return tmp<Field<Type> >
        (
            new Field<Type>
            (
                ListListOps::combine<Field<Type> >
                (
                    allValues,
                    accessOp<Field<Type> >()
                )
            )
        );
    }
    else
    {
        return field;
    }
}


// ************************************************************************* //
