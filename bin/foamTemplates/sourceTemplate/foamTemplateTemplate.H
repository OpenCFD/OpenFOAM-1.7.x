/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::ClassName

Description

SourceFiles
    ClassNameI.H
    ClassName.C
    ClassNameIO.C

\*---------------------------------------------------------------------------*/

#ifndef ClassName_H
#define ClassName_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class someClass;


// Forward declaration of friend functions and operators

template<TemplateClassArgument>
class ClassName;

template<TemplateClassArgument>
Istream& operator>>(Istream&, ClassName<TemplateArgument>&);

template<TemplateClassArgument>
Ostream& operator<<(Ostream&, const ClassName<TemplateArgument>&);


/*---------------------------------------------------------------------------*\
                         Class ClassName Declaration
\*---------------------------------------------------------------------------*/

template<TemplateClassArgument>
class ClassName
:
    public baseClassName
{
    // Private data

        dataType data_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        ClassName(const ClassName<TemplateArgument>&);

        //- Disallow default bitwise assignment
        void operator=(const ClassName<TemplateArgument>&);


public:

    // Static data members

        //- Static data someStaticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        ClassName();

        //- Construct from components
        ClassName(const dataType& data);

        //- Construct from Istream
        ClassName(Istream&);

        //- Construct as copy
        ClassName(const ClassName<TemplateArgument>&);


    // Selectors

        //- Select null constructed
        static autoPtr<ClassName<TemplateArgument> > New();


    //- Destructor
    ~ClassName();


    // Member Functions

    // Member Operators

        void operator=(const ClassName<TemplateArgument>&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>> <TemplateArgument>
        (Istream&, ClassName<TemplateArgument>&);

        friend Ostream& operator<< <TemplateArgument>
        (Ostream&, const ClassName<TemplateArgument>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ClassNameI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ClassName.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
