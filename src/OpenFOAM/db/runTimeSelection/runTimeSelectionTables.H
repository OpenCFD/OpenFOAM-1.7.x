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
    Foam::runTimeSelectionTables

Description
    Macros to enable the easy declaration of run-time selection tables.

    declareRunTimeSelectionTable is used to create a run-time selection table
    for a base-class which holds constructor pointers on the table.

    declareRunTimeNewSelectionTable is used to create a run-time selection
    table for a derived-class which holds "New" pointers on the table.

\*---------------------------------------------------------------------------*/

#include "token.H"

#ifndef runTimeSelectionTables_H
#define runTimeSelectionTables_H

#include "autoPtr.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// external use:
// ~~~~~~~~~~~~~
// declare a run-time selection:
#define declareRunTimeSelectionTable\
(autoPtr,baseType,argNames,argList,parList)                                   \
                                                                              \
    /* Construct from argList function pointer type */                        \
    typedef autoPtr< baseType > (*argNames##ConstructorPtr)argList;           \
                                                                              \
    /* Construct from argList function table type */                          \
    typedef HashTable< argNames##ConstructorPtr, word, string::hash >         \
        argNames##ConstructorTable;                                           \
                                                                              \
    /* Construct from argList function pointer table pointer */               \
    static argNames##ConstructorTable* argNames##ConstructorTablePtr_;        \
                                                                              \
    /* Class to add constructor from argList to table */                      \
    template< class baseType##Type >                                          \
    class add##argNames##ConstructorToTable                                   \
    {                                                                         \
    public:                                                                   \
                                                                              \
        static autoPtr< baseType > New argList                                \
        {                                                                     \
            return autoPtr< baseType >(new baseType##Type parList);           \
        }                                                                     \
                                                                              \
        add##argNames##ConstructorToTable                                     \
        (                                                                     \
            const word& lookup = baseType##Type::typeName                     \
        )                                                                     \
        {                                                                     \
            construct##argNames##ConstructorTables();                         \
            argNames##ConstructorTablePtr_->insert(lookup, New);              \
        }                                                                     \
                                                                              \
        ~add##argNames##ConstructorToTable()                                  \
        {                                                                     \
            destroy##argNames##ConstructorTables();                           \
        }                                                                     \
    };                                                                        \
                                                                              \
    /* Table constructor called from the table add function */                \
    static void construct##argNames##ConstructorTables();                     \
                                                                              \
    /* Table destructor called from the table add function destructor */      \
    static void destroy##argNames##ConstructorTables()


// external use:
// ~~~~~~~~~~~~~
// declare a run-time selection for derived classes:
#define declareRunTimeNewSelectionTable\
(autoPtr,baseType,argNames,argList,parList)                                   \
                                                                              \
    /* Construct from argList function pointer type */                        \
    typedef autoPtr< baseType > (*argNames##ConstructorPtr)argList;           \
                                                                              \
    /* Construct from argList function table type */                          \
    typedef HashTable< argNames##ConstructorPtr, word, string::hash >         \
        argNames##ConstructorTable;                                           \
                                                                              \
    /* Construct from argList function pointer table pointer */               \
    static argNames##ConstructorTable* argNames##ConstructorTablePtr_;        \
                                                                              \
    /* Class to add constructor from argList to table */                      \
    template< class baseType##Type >                                          \
    class add##argNames##ConstructorToTable                                   \
    {                                                                         \
    public:                                                                   \
                                                                              \
        static autoPtr< baseType > New##baseType argList                      \
        {                                                                     \
            return autoPtr< baseType >(baseType##Type::New parList.ptr());    \
        }                                                                     \
                                                                              \
        add##argNames##ConstructorToTable                                     \
        (                                                                     \
            const word& lookup = baseType##Type::typeName                     \
        )                                                                     \
        {                                                                     \
            construct##argNames##ConstructorTables();                         \
            argNames##ConstructorTablePtr_->insert                            \
            (                                                                 \
                lookup,                                                       \
                New##baseType                                                 \
            );                                                                \
        }                                                                     \
                                                                              \
        ~add##argNames##ConstructorToTable()                                  \
        {                                                                     \
            destroy##argNames##ConstructorTables();                           \
        }                                                                     \
    };                                                                        \
                                                                              \
    /* Table constructor called from the table add function */                \
    static void construct##argNames##ConstructorTables();                     \
                                                                              \
    /* Table destructor called from the table add function destructor */      \
    static void destroy##argNames##ConstructorTables()


// internal use:
// constructor aid
#define defineRunTimeSelectionTableConstructor\
(baseType,argNames)                                                           \
                                                                              \
    /* Table constructor called from the table add function */                \
    void baseType::construct##argNames##ConstructorTables()                   \
    {                                                                         \
        static bool constructed = false;                                      \
                                                                              \
        if (!constructed)                                                     \
        {                                                                     \
            baseType::argNames##ConstructorTablePtr_                          \
                = new baseType::argNames##ConstructorTable;                   \
                                                                              \
            constructed = true;                                               \
        }                                                                     \
    }


// internal use:
// destructor aid
#define defineRunTimeSelectionTableDestructor\
(baseType,argNames)                                                           \
                                                                              \
    /* Table destructor called from the table add function destructor */      \
    void baseType::destroy##argNames##ConstructorTables()                     \
    {                                                                         \
        if (baseType::argNames##ConstructorTablePtr_)                         \
        {                                                                     \
            delete baseType::argNames##ConstructorTablePtr_;                  \
            baseType::argNames##ConstructorTablePtr_ = NULL;                  \
        }                                                                     \
    }


// internal use:
// create pointer to hash-table of functions
#define defineRunTimeSelectionTablePtr\
(baseType,argNames)                                                           \
                                                                              \
    /* Define the constructor function table */                               \
    baseType::argNames##ConstructorTable*                                     \
        baseType::argNames##ConstructorTablePtr_ = NULL


// not much in use:
#define defineTemplateRunTimeSelectionTablePtr(baseType,argNames)             \
                                                                              \
    /* Define the constructor function table */                               \
    typename baseType::argNames##ConstructorTable*                            \
        baseType::argNames##ConstructorTablePtr_ = NULL


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// external use:
// ~~~~~~~~~~~~~
// define run-time selection table
#define defineRunTimeSelectionTable\
(baseType,argNames)                                                           \
                                                                              \
    defineRunTimeSelectionTablePtr(baseType,argNames);                        \
    defineRunTimeSelectionTableConstructor(baseType,argNames);                \
    defineRunTimeSelectionTableDestructor(baseType,argNames)


// external use:
// ~~~~~~~~~~~~~
// define run-time selection table for template classes
// use when baseType doesn't need a template argument (eg, is a typedef)
#define defineTemplateRunTimeSelectionTable\
(baseType,argNames)                                                           \
                                                                              \
    template<>                                                                \
    defineRunTimeSelectionTablePtr(baseType,argNames);                        \
    template<>                                                                \
    defineRunTimeSelectionTableConstructor(baseType,argNames);                \
    template<>                                                                \
    defineRunTimeSelectionTableDestructor(baseType,argNames)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// internal use:
// constructor aid
// use when baseType requires the Targ template argument
#define defineTemplatedRunTimeSelectionTableConstructor\
(baseType,argNames,Targ)                                                      \
                                                                              \
    /* Table constructor called from the table add function */                \
    void baseType< Targ >::construct##argNames##ConstructorTables()           \
    {                                                                         \
        static bool constructed = false;                                      \
                                                                              \
        if (!constructed)                                                     \
        {                                                                     \
            baseType< Targ >::argNames##ConstructorTablePtr_                  \
                = new baseType< Targ >::argNames##ConstructorTable;           \
                                                                              \
            constructed = true;                                               \
        }                                                                     \
    }


// internal use:
// destructor aid
// use when baseType requires the Targ template argument
#define defineTemplatedRunTimeSelectionTableDestructor\
(baseType,argNames,Targ)                                                      \
                                                                              \
    /* Table destructor called from the table add function destructor */      \
    void baseType< Targ >::destroy##argNames##ConstructorTables()             \
    {                                                                         \
        if (baseType< Targ >::argNames##ConstructorTablePtr_)                 \
        {                                                                     \
            delete baseType< Targ >::argNames##ConstructorTablePtr_;          \
            baseType< Targ >::argNames##ConstructorTablePtr_ = NULL;          \
        }                                                                     \
    }


// internal use:
// create pointer to hash-table of functions
// use when baseType requires the Targ template argument
#define defineTemplatedRunTimeSelectionTablePtr\
(baseType,argNames,Targ)                                                      \
                                                                              \
    /* Define the constructor function table */                               \
    baseType< Targ >::argNames##ConstructorTable*                             \
        baseType< Targ >::argNames##ConstructorTablePtr_ = NULL


// external use:
// ~~~~~~~~~~~~~
// define run-time selection table for template classes
// use when baseType requires the Targ template argument
#define defineTemplatedRunTimeSelectionTable\
(baseType,argNames,Targ)                                                      \
                                                                              \
    template<>                                                                \
    defineTemplatedRunTimeSelectionTablePtr(baseType,argNames,Targ);          \
    template<>                                                                \
    defineTemplatedRunTimeSelectionTableConstructor(baseType,argNames,Targ);  \
    template<>                                                                \
    defineTemplatedRunTimeSelectionTableDestructor(baseType,argNames,Targ)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
