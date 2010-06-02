/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Matcher, class StringType>
Foam::labelList Foam::findMatchingStrings
(
    const Matcher& matcher,
    const UList<StringType>& lst,
    const bool invert
)
{
    labelList indices(lst.size());

    label nElem = 0;
    forAll(lst, elemI)
    {
        if (matcher.match(lst[elemI]) ? !invert : invert)
        {
            indices[nElem++] = elemI;
        }
    }
    indices.setSize(nElem);

    return indices;
}


template<class Matcher, class StringListType>
StringListType Foam::subsetMatchingStrings
(
    const Matcher& matcher,
    const StringListType& lst,
    const bool invert
)
{
    StringListType newLst(lst.size());

    label nElem = 0;
    forAll(lst, elemI)
    {
        if (matcher.match(lst[elemI]) ? !invert : invert)
        {
            newLst[nElem++] = lst[elemI];
        }
    }
    newLst.setSize(nElem);

    return newLst;
}


template<class Matcher, class StringListType>
void Foam::inplaceSubsetMatchingStrings
(
    const Matcher& matcher,
    StringListType& lst,
    const bool invert
)
{
    label nElem = 0;
    forAll(lst, elemI)
    {
        if (matcher.match(lst[elemI]) ? !invert : invert)
        {
            lst[nElem++] = lst[elemI];
        }
    }
    lst.setSize(nElem);
}


// ************************************************************************* //
