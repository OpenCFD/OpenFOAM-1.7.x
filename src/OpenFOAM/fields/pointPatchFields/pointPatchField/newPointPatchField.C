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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type> > Foam::pointPatchField<Type>::New
(
    const word& patchFieldType,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
{
    if (debug)
    {
        Info<< "PointPatchField<Type>::"
               "New(const word&, const pointPatch&, const Field<Type>&) : "
               "constructing pointPatchField<Type>"
            << endl;
    }

    typename pointPatchConstructorTable::iterator cstrIter =
        pointPatchConstructorTablePtr_->find(patchFieldType);

    if (cstrIter == pointPatchConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "PointPatchField<Type>::New"
            "(const word&, const pointPatch&, const Field<Type>&)"
        )   << "Unknown patchTypefield type "
            << patchFieldType
            << endl << endl
            << "Valid patchField types are :" << endl
            << pointPatchConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    typename pointPatchConstructorTable::iterator patchTypeCstrIter =
        pointPatchConstructorTablePtr_->find(p.type());

    if (patchTypeCstrIter != pointPatchConstructorTablePtr_->end())
    {
        return autoPtr<pointPatchField<Type> >(patchTypeCstrIter()(p, iF));
    }
    else
    {
        return autoPtr<pointPatchField<Type> >(cstrIter()(p, iF));
    }
}


template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type> > Foam::pointPatchField<Type>::New
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
{
    if (debug)
    {
        Info<< "PointPatchField<Type>::"
               "New(const pointPatch&, const Field<Type>&, const dictionary&)"
               " : constructing pointPatchField<Type>"
            << endl;
    }

    word patchFieldType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter
        = dictionaryConstructorTablePtr_->find(patchFieldType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        if (!disallowGenericPointPatchField)
        {
            cstrIter = dictionaryConstructorTablePtr_->find("generic");
        }

        if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalIOErrorIn
            (
                "PointPatchField<Type>::"
                "New(const pointPatch&, const Field<Type>&, const dictionary&)",
                dict
            )   << "Unknown patchField type " << patchFieldType
                << " for patch type " << p.type() << endl << endl
                << "Valid patchField types are :" << endl
                << dictionaryConstructorTablePtr_->toc()
                << exit(FatalIOError);
        }
    }

    if
    (
       !dict.found("patchType")
     || word(dict.lookup("patchType")) != p.type()
    )
    {
        typename dictionaryConstructorTable::iterator patchTypeCstrIter
            = dictionaryConstructorTablePtr_->find(p.type());

        if
        (
            patchTypeCstrIter != dictionaryConstructorTablePtr_->end()
         && patchTypeCstrIter() != cstrIter()
        )
        {
            FatalIOErrorIn
            (
                "PointPatchField<Type>const pointPatch&, "
                "const Field<Type>&, const dictionary&)",
                dict
            )   << "inconsistent patch and patchField types for \n"
                << "    patch type " << p.type()
                << " and patchField type " << patchFieldType
                << exit(FatalIOError);
        }
    }

    return autoPtr<pointPatchField<Type> >(cstrIter()(p, iF, dict));
}


// Return a pointer to a new patch created on freestore from
// a given pointPatchField<Type> mapped onto a new patch
template<class Type>
Foam::autoPtr<Foam::pointPatchField<Type> > Foam::pointPatchField<Type>::New
(
    const pointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& pfMapper
)
{
    if (debug)
    {
        Info<< "PointPatchField<Type>::"
               "New(const pointPatchField<Type>&,"
               " const pointPatch&, const Field<Type>&, "
               "const pointPatchFieldMapper&) : "
               "constructing pointPatchField<Type>"
            << endl;
    }

    typename patchMapperConstructorTable::iterator cstrIter =
        patchMapperConstructorTablePtr_->find(ptf.type());

    if (cstrIter == patchMapperConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "PointPatchField<Type>::"
            "New(const pointPatchField<Type>&, "
            "const pointPatch&, const Field<Type>&, "
            "const pointPatchFieldMapper&)"
        )   << "unknown patchTypefield type "
            << ptf.type() << endl << endl
            << "Valid patchField types are :" << endl
            << patchMapperConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<pointPatchField<Type> >(cstrIter()(ptf, p, iF, pfMapper));
}


// ************************************************************************* //
