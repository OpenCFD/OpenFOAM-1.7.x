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

InNamspace
    Foam

Description
    Operations on lists of strings.

SourceFiles
    stringListOpsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef stringListOps_H
#define stringListOps_H

#include "labelList.H"
#include "stringList.H"
#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //- Return list indices for strings matching the regular expression
    template<class StringType>
    labelList findStrings
    (
        const char* regexpPattern,
        const UList<StringType>&
    );

    //- Return list indices for strings matching the regular expression
    template<class StringType>
    labelList findStrings
    (
        const std::string& regexpPattern,
        const UList<StringType>&
    );

    //- Return list indices for strings matching the regular expression
    template<class StringType>
    labelList findStrings
    (
        const wordRe&,
        const UList<StringType>&
    );

    //- Return list indices for strings matching one of the regular expression
    template<class StringType>
    labelList findStrings
    (
        const UList<wordRe>&,
        const UList<StringType>&
    );

    //- Return true if string matches one of the regular expressions
    template<class StringType>
    bool findStrings
    (
        const UList<wordRe>&,
        const StringType& str
    );

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "stringListOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
