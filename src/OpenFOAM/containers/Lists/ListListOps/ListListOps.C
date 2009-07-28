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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template <class AccessType, class T, class AccessOp>
AccessType ListListOps::combine(const List<T>& lst, AccessOp aop)
{
    label sum = 0;

    forAll(lst, lstI)
    {
        sum += aop(lst[lstI]).size();
    }

    AccessType result(sum);

    label globalElemI = 0;

    forAll(lst, lstI)
    {
        const T& sub = lst[lstI];

        forAll(aop(sub), elemI)
        {
            result[globalElemI++] = aop(sub)[elemI];
        }
    }
    return result;
}


template <class T, class AccessOp>
labelList ListListOps::subSizes(const List<T>& lst, AccessOp aop)
{
    labelList sizes(lst.size());

    forAll(lst, lstI)
    {
        sizes[lstI] = aop(lst[lstI]).size();
    }
    return sizes;
}


template <class AccessType, class T, class AccessOp, class OffsetOp>
AccessType ListListOps::combineOffset
(
    const List<T>& lst,
    const labelList& sizes,
    AccessOp aop,
    OffsetOp oop
)
{
    label sum = 0;

    forAll(lst, lstI)
    {
        sum += aop(lst[lstI]).size();
    }

    AccessType result(sum);

    label globalElemI = 0;

    label offset = 0;

    forAll(lst, lstI)
    {
        const T& sub = lst[lstI];

        forAll(aop(sub), elemI)
        {
            result[globalElemI++] = oop(aop(sub)[elemI], offset);
        }

        offset += sizes[lstI];
    }
    return result;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
