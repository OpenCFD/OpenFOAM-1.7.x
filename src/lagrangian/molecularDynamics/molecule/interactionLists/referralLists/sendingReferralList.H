/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    Foam::sendingReferralList

Description

SourceFiles
    sendingReferralListI.H
    sendingReferralList.C
    sendingReferralListIO.C

\*---------------------------------------------------------------------------*/

#ifndef sendingReferralList_H
#define sendingReferralList_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class sendingReferralList Declaration
\*---------------------------------------------------------------------------*/

class sendingReferralList
:
    public labelList
{
    // Private data

        label destinationProc_;


public:

    // Constructors

        //- Construct null
        sendingReferralList();

        //- Construct from components
        sendingReferralList
        (
            const label destinationProc,
            const labelList& cellsToSend
        );

        //- Construct as copy
        sendingReferralList(const sendingReferralList&);


    // Destructor

        ~sendingReferralList();


    // Member Functions

        // Access

            inline label destinationProc() const;


    // Member Operators

        void operator=(const sendingReferralList&);


    // Friend Operators

        friend bool operator==
        (
            const sendingReferralList& a,
            const sendingReferralList& b
        );

        inline friend bool operator!=
        (
            const sendingReferralList& a,
            const sendingReferralList& b
        );


    // IOstream Operators

        friend Istream& operator>>
        (
            Istream&,
            sendingReferralList&
        );

        friend Ostream& operator<<
        (
            Ostream&,
            const sendingReferralList&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sendingReferralListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
