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

\*---------------------------------------------------------------------------*/

#include "surfaceWriter.H"

#include "MeshedSurfaceProxy.H"
#include "nullSurfaceWriter.H"
#include "proxySurfaceWriter.H"

#include "HashTable.H"
#include "word.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr< Foam::surfaceWriter<Type> >
Foam::surfaceWriter<Type>::New(const word& writeType)
{
    typename wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(writeType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        // not supported for this data type, but it generally does work
        // (it handles the 'bool' specialization - ie, geometry write)
        if
        (
            Foam::surfaceWriter<bool>::wordConstructorTablePtr_->found
            (
                writeType
            )
        )
        {
            // use 'null' handler instead
            return autoPtr< surfaceWriter<Type> >
            (
                new nullSurfaceWriter<Type>()
            );
        }
        else if (MeshedSurfaceProxy<face>::canWriteType(writeType))
        {
            // generally unknown, but can be written via MeshedSurfaceProxy
            // use 'proxy' handler instead
            return autoPtr< surfaceWriter<Type> >
            (
                new proxySurfaceWriter<Type>(writeType)
            );
        }

        if (cstrIter == wordConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "surfaceWriter::New(const word&)"
            )   << "Unknown write type \"" << writeType << "\"\n\n"
                << "Valid write types : "
                << wordConstructorTablePtr_->toc() << nl
                << "Valid proxy types : "
                << MeshedSurfaceProxy<face>::writeTypes() << endl
                << exit(FatalError);
        }
    }

    return autoPtr< surfaceWriter<Type> >(cstrIter()());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::surfaceWriter<Type>::surfaceWriter()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::surfaceWriter<Type>::~surfaceWriter()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// ************************************************************************* //
