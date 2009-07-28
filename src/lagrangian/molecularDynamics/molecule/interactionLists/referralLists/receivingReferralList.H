/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    Foam::receivingReferralList

Description

SourceFiles
    receivingReferralListI.H
    receivingReferralList.C
    receivingReferralListIO.C

\*---------------------------------------------------------------------------*/

#ifndef receivingReferralList_H
#define receivingReferralList_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class receivingReferralList Declaration
\*---------------------------------------------------------------------------*/

class receivingReferralList
:
    public labelListList
{
    // Private data

        label sourceProc_;


public:

    // Constructors

        //- Construct null
        receivingReferralList();

        //- Construct from components
        receivingReferralList
        (
            const label sourceProc,
            const labelListList& refCellsToSendTo
        );

        //- Construct as copy
        receivingReferralList(const receivingReferralList&);


    // Destructor

        ~receivingReferralList();


    // Member Functions

        // Access

            inline label sourceProc() const;


    // Member Operators

        void operator=(const receivingReferralList&);


    // Friend Operators

        friend bool operator==
        (
            const receivingReferralList& a,
            const receivingReferralList& b
        );

        inline friend bool operator!=
        (
            const receivingReferralList& a,
            const receivingReferralList& b
        );


    // IOstream Operators

        friend Istream& operator>>
        (
            Istream&,
            receivingReferralList&
        );

        friend Ostream& operator<<
        (
            Ostream&,
            const receivingReferralList&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "receivingReferralListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
