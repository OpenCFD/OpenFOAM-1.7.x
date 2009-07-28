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
    Foam::chemistryReader

Description
    Abstract class for reading chemistry

SourceFiles
    chemistryReader.C

\*---------------------------------------------------------------------------*/

#ifndef chemistryReader_H
#define chemistryReader_H

#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "Reaction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class chemistryReader Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class chemistryReader
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        chemistryReader(const chemistryReader&);

        //- Disallow default bitwise assignment
        void operator=(const chemistryReader&);


public:

    //- Runtime type information
    TypeName("chemistryReader");

    //- The type of thermo package the reader was instantiated for
    typedef ThermoType thermoType;


    // Constructors

        //- Construct null
        chemistryReader()
        {}


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            chemistryReader,
            dictionary,
            (
                const dictionary& thermoDict
            ),
            (thermoDict)
        );


    // Selectors

        //- Select constructed from dictionary
        static autoPtr<chemistryReader> New(const dictionary& thermoDict);


    //- Destructor
    virtual ~chemistryReader()
    {}


    // Member Functions

        virtual const speciesTable& species() const = 0;

        virtual const HashPtrTable<ThermoType>& speciesThermo() const = 0;

        virtual const SLPtrList<Reaction<ThermoType> >& reactions() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChemistryReader(Thermo)                                           \
                                                                              \
defineTemplateTypeNameAndDebug                                                \
(                                                                             \
    chemistryReader<Thermo>,                                                  \
    0                                                                         \
);                                                                            \
                                                                              \
defineTemplateRunTimeSelectionTable(chemistryReader<Thermo>, dictionary);


#define makeChemistryReaderType(SS, Thermo)                                   \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<Thermo>, 0);                       \
                                                                              \
    chemistryReader<Thermo>::adddictionaryConstructorToTable<SS<Thermo> >     \
        add##SS##Thermo##ConstructorToTable_;


#define addChemistryReaderType(SS, Thermo)                                    \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS, 0);                               \
                                                                              \
    chemistryReader<Thermo>::adddictionaryConstructorToTable<SS>              \
        add##SS##Thermo##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "chemistryReader.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
