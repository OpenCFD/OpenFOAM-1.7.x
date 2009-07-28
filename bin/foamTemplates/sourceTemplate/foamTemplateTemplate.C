/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "ClassName.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<TemplateClassArgument>
const dataType Foam::ClassName<TemplateArgument>::staticData();


// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<TemplateClassArgument>
Foam::ClassName<TemplateArgument>::ClassName()
:
    baseClassName(),
    data_()
{}


template<TemplateClassArgument>
Foam::ClassName<TemplateArgument>::ClassName(const dataType& data)
:
    baseClassName(),
    data_(data)
{}


template<TemplateClassArgument>
Foam::ClassName<TemplateArgument>::ClassName
(
    const ClassName<TemplateArgument>&
)
:
    baseClassName(),
    data_()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<TemplateClassArgument>
Foam::autoPtr<Foam::ClassName<TemplateArgument> >
Foam::ClassName<TemplateArgument>::New()
{
    return autoPtr<ClassName<TemplateArgument> >
    (
        new ClassName<TemplateArgument>
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<TemplateClassArgument>
Foam::ClassName<TemplateArgument>::~ClassName()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<TemplateClassArgument>
void Foam::ClassName<TemplateArgument>::operator=
(
    const ClassName<TemplateArgument>& rhs
)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::ClassName<TemplateArgument>::operator="
            "(const Foam::ClassName<TemplateArgument>&)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
