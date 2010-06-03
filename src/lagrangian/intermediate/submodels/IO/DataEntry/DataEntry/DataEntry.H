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

Class
    Foam::DataEntry

Description
    Top level data entry class for use in dictionaries. Provides a mechanism
    to specify a variable as a certain type, e.g. constant or table, and
    provide functions to return the (interpolated) value, and integral between
    limits.

SourceFiles
    DataEntry.C
    NewDataEntry.C

\*---------------------------------------------------------------------------*/

#ifndef DataEntry_H
#define DataEntry_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class DataEntry;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const DataEntry<Type>&
);

/*---------------------------------------------------------------------------*\
                         Class DataEntry Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DataEntry
{
    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const DataEntry<Type>&);


protected:

    // Protected data

        //- Name of entry
        const word name_;


public:

    //- Runtime type information
    TypeName("DataEntry")

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        DataEntry,
        dictionary,
        (
            const word& entryName,
            Istream& is
        ),
        (entryName, is)
    );


    // Constructor

        //- Construct from entry name
        DataEntry(const word& entryName);

        //- Copy constructor
        DataEntry(const DataEntry<Type>& de);


    //- Selector
    static autoPtr<DataEntry<Type> > New
    (
        const word& entryName,
        const dictionary& dict
    );


    //- Destructor
    virtual ~DataEntry();


    // Member Functions

        // Access

            //- Return the name of the entry
            const word& name() const;


        // Evaluation

            //- Return value as a function of (scalar) independent variable
            virtual Type value(const scalar x) const = 0;

            //- Integrate between two (scalar) values
            virtual Type integrate(const scalar x1, const scalar x2) const = 0;


    //- Ostream Operator
    friend Ostream& operator<< <Type>
    (
        Ostream&,
        const DataEntry<Type>&
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDataEntry(Type)                                                   \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(DataEntry<Type>, 0);                  \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        DataEntry<Type>,                                                      \
        dictionary                                                            \
    );


#define makeDataEntryType(SS, Type)                                           \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                         \
                                                                              \
    DataEntry<Type>::adddictionaryConstructorToTable<SS<Type> >               \
        add##SS##Type##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DataEntry.C"
#   include "NewDataEntry.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
