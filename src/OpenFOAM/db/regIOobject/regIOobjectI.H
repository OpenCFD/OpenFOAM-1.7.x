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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::regIOobject::ownedByRegistry() const
{
    return ownedByRegistry_;
}


inline void Foam::regIOobject::store()
{
    ownedByRegistry_ = true;
}


template<class Type>
inline Type& Foam::regIOobject::store(Type* tPtr)
{
    if (!tPtr)
    {
        FatalErrorIn("Type& regIOobject::store(Type*)")
            << "object deallocated"
            << abort(FatalError);
    }

    tPtr->regIOobject::ownedByRegistry_ = true;

    return *tPtr;
}


template<class Type>
inline Type& Foam::regIOobject::store(autoPtr<Type>& atPtr)
{
    Type* tPtr = atPtr.ptr();

    if (!tPtr)
    {
        FatalErrorIn
            (
                "Type& regIOobject::store(autoPtr<Type>&)"
            )   << "object deallocated"
                << abort(FatalError);
    }

    tPtr->regIOobject::ownedByRegistry_ = true;

    return *tPtr;
}


inline void Foam::regIOobject::release()
{
    ownedByRegistry_ = false;
}


inline Foam::label Foam::regIOobject::eventNo() const
{
    return eventNo_;
}

inline Foam::label& Foam::regIOobject::eventNo()
{
    return eventNo_;
}


// ************************************************************************* //
