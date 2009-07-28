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
    Foam::IntegrationScheme

Description
    Top level model for Integration schemes

SourceFiles
    IntegrationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef IntegrationScheme_H
#define IntegrationScheme_H

#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class IntegrationScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class IntegrationScheme
{

public:

    //- Helper class to supply results of integration
    class integrationResult
    {
        //- Integration value
        Type value_;

        //- Average value across integration step
        Type average_;


    public:

        //- Constructor
        integrationResult()
        :
            value_(pTraits<Type>::zero),
            average_(pTraits<Type>::zero)
        {}


        // Member functions

            // Access

                //- Return const access to the value
                Type value() const
                {
                    return value_;
                }

                //- Return const access to the average
                Type average() const
                {
                    return average_;
                }


            // Edit

                //- Return access to the value for changing
                Type& value()
                {
                    return value_;
                }

                //- Return access to the average for changing
                Type& average()
                {
                    return average_;
                }
        };


private:

    // Private data

        //- Name of the Integration variable
        const word& phiName_;

        //- Reference to the dictionary
        const dictionary& dict_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        IntegrationScheme(const IntegrationScheme&);

        //- Disallow default bitwise assignment
        void operator=(const IntegrationScheme&);


public:

    //- Runtime type information
    TypeName("IntegrationScheme");


    //- Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             IntegrationScheme,
             dictionary,
             (
                 const word& phiName,
                 const dictionary& dict
             ),
             (phiName, dict)
         );


    // Constructors

        //- Construct from components
        IntegrationScheme(const word& phiName, const dictionary& dict);


    // Selectors

        //- Return a reference to the selected radiation model
        static autoPtr<IntegrationScheme> New
        (
            const word& phiName,
            const dictionary& dict
        );


    //- Destructor
    virtual ~IntegrationScheme();


    // Member Functions

        //- Perform the Integration
        virtual integrationResult integrate
        (
            const Type phi,
            const scalar dt,
            const Type alpha,
            const scalar beta
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeIntegrationScheme(Type)                                           \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(IntegrationScheme<Type>, 0);          \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        IntegrationScheme<Type>,                                               \
        dictionary                                                            \
    );


#define makeIntegrationSchemeType(SS, Type)                                   \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                         \
                                                                              \
    IntegrationScheme<Type>::adddictionaryConstructorToTable<SS<Type> >       \
        add##SS##Type##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IntegrationScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

