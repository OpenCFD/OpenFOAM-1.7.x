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

Class
    Foam::pointSourceProperties

Description
    Helper class to describe point source properties

SourceFiles
    pointSourceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef pointSourceProperties_H
#define pointSourceProperties_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "Time.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class pointSourceProperties Declaration
\*---------------------------------------------------------------------------*/

class pointSourceProperties
{
protected:

    // Protected data

        typedef Tuple2<word, scalar> fieldNameValuePair;

        //- Source name
        word name_;

        //- Time start
        scalar timeStart_;

        //- Duration
        scalar duration_;

        //- Point location
        point location_;

        //- List of source field name vs value pairs
        List<fieldNameValuePair> fieldData_;


public:

    // Constructors

        //- Construct null
        pointSourceProperties();

        //- Construct from dictionary
        pointSourceProperties(const dictionary& dict);

        //- Construct from Istream
        pointSourceProperties(Istream& is);

        //- Copy constructor
        pointSourceProperties(const pointSourceProperties&);


    // Member Functions

        // Access

            //- Return const access to the source name
            inline const word& name() const;

            //- Return const access to the time start
            inline scalar timeStart() const;

            //- Return const access to the time end
            inline scalar timeEnd() const;

            //- Return const access to the duration
            inline scalar duration() const;

            //- Return const access to the point location
            inline const point& location() const;

            //- Return const access to the source field name vs value pairs
            inline const List<fieldNameValuePair>& fieldData() const;


        // Edit

            //- Return access to the source name
            inline word& name();

            //- Return access to the time start
            inline scalar& timeStart();

            //- Return access to the duration
            inline scalar& duration();

            //- Return access to the point location
            inline point& location();

            //- Return access to the source field name vs value pairs
            inline List<fieldNameValuePair>& fieldData();


        // Member Operators

            void operator=(const pointSourceProperties&);

        // IOstream operators

            friend Istream& operator>>(Istream&, pointSourceProperties&);
            friend Ostream& operator<<(Ostream&, const pointSourceProperties&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointSourcePropertiesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
