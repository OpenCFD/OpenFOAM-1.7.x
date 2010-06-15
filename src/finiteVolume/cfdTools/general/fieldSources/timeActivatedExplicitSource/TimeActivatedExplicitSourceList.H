/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    Foam::TimeActivatedExplicitSourceList

Description
    List of time activeted explict sources

SourceFiles
    TimeActivatedExplicitSourceList.C

\*---------------------------------------------------------------------------*/

#ifndef TimeActivatedExplicitSourceList_H
#define TimeActivatedExplicitSourceList_H

#include "IOPtrList.H"
#include "TimeActivatedExplicitSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

template<class Type>
class TimeActivatedExplicitSource;

template<class Type>
class TimeActivatedExplicitSourceList;

// Forward declaration of friend functions

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const TimeActivatedExplicitSourceList<Type>&
);

/*---------------------------------------------------------------------------*\
             Class TimeActivatedExplicitSourceList Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class TimeActivatedExplicitSourceList
:
    public IOPtrList<TimeActivatedExplicitSource<Type> >
{
private:

    // Private data

        //- Name of source properties list
        word name_;

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Dimensions of source properties
        dimensionSet dimensions_;

        //- List of field names the source is working on
        wordList fieldNames_;


    // Private member functions

        //- Disallow default bitwise copy construct
        TimeActivatedExplicitSourceList
        (
            const TimeActivatedExplicitSourceList<Type>&
        );

        //- Disallow default bitwise assignment
        void operator=(const TimeActivatedExplicitSourceList<Type>&);


public:

    // Constructors

        //- Construct from components with list of field names
        TimeActivatedExplicitSourceList
        (
            const word& name,
            const fvMesh& mesh,
            const dimensionSet& dimensions,
            const wordList& fieldNames
        );

        //- Construct from components with single field name
        TimeActivatedExplicitSourceList
        (
            const word& name,
            const fvMesh& mesh,
            const dimensionSet& dimensions,
            const word& fieldName
        );


    // Member Functions

        // Evaluation

            //- Return the source for field, fieldI
            tmp<DimensionedField<Type, volMesh> > Su
            (
                const label fieldI = 0
            );

            //- Return the total source for all fields
            tmp<DimensionedField<Type, volMesh> > SuTot();



        // I-O

            //- Read data from Istream
            bool readData(Istream& is);

            //- Write data to Istream
            bool writeData(Ostream& os) const;

            //- Ostream operator
            friend Ostream& operator<< <Type>
            (
                Ostream& os,
                const TimeActivatedExplicitSourceList& sources
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "TimeActivatedExplicitSourceList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
