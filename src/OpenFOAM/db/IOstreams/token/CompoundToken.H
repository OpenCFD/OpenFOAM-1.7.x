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

Class
    Foam::CompoundToken

Description
    An abstract base class for managing compound tokens

\*---------------------------------------------------------------------------*/

#ifndef CompoundToken_H
#define CompoundToken_H

#include "refCount.H"
#include "typeInfo.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class CompoundToken;
Ostream& operator<<(Ostream&, const CompoundToken&);


/*---------------------------------------------------------------------------*\
                           Class CompoundToken Declaration
\*---------------------------------------------------------------------------*/

class CompoundToken
:
    public refCount
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        CompoundToken(const CompoundToken&);

        //- Disallow default bitwise assignment
        void operator=(const CompoundToken&);


public:

    //- Runtime type information
    virtual const word& type() const = 0;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            CompoundToken,
            Istream,
            (const word& type, Istream& is),
            (type, is)
        );


    // Constructors

        //- Construct null
        CompoundToken()
        {}

        //- Return the clone as this and increment reference count
        virtual autoPtr<CompoundToken> clone() const = 0;


    // Selectors

        //- Select null constructed
        static autoPtr<CompoundToken> New(const word& type, Istream& is);


    // Destructor

        virtual ~CompoundToken();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write

            virtual void write(Istream&) = 0;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const CompoundToken&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
