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

#include "genericPointPatchField.H"
#include "pointPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
genericPointPatchField<Type>::genericPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    calculatedPointPatchField<Type>(p, iF)
{
    notImplemented
    (
        "genericPointPatchField<Type>::genericPointPatchField"
        "(const pointPatch& p, const DimensionedField<Type, volMesh>& iF)"
    );
}


template<class Type>
genericPointPatchField<Type>::genericPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    calculatedPointPatchField<Type>(p, iF, dict),
    actualTypeName_(dict.lookup("type")),
    dict_(dict)
{
    for
    (
        dictionary::const_iterator iter = dict_.begin();
        iter != dict_.end();
        ++iter
    )
    {
        if (iter().keyword() != "type")
        {
            if
            (
                iter().isStream()
             && iter().stream().size()
            )
            {
                ITstream& is = iter().stream();

                // Read first token
                token firstToken(is);

                if
                (
                    firstToken.isWord()
                 && firstToken.wordToken() == "nonuniform"
                )
                {
                    token fieldToken(is);

                    if (!fieldToken.isCompound())
                    {
                        if
                        (
                            fieldToken.isLabel()
                         && fieldToken.labelToken() == 0
                        )
                        {
                            scalarFields_.insert
                            (
                                iter().keyword(),
                                new scalarField(0)
                            );
                        }
                        else
                        {
                            FatalIOErrorIn
                            (
                                "genericPointPatchField<Type>::"
                                "genericPointPatchField"
                                "(const pointPatch&, const Field<Type>&, "
                                "const dictionary&)",
                                dict
                            )   << "\n    token following 'nonuniform' "
                                  "is not a compound"
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->dimensionedInternalField().name()
                                << " in file "
                                << this->dimensionedInternalField().objectPath()
                            << exit(FatalIOError);
                        }
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<scalar> >::typeName
                    )
                    {
                        scalarField* fPtr = new scalarField;
                        fPtr->transfer
                        (
                            dynamicCast<token::Compound<List<scalar> > >
                            (
                                fieldToken.transferCompoundToken()
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorIn
                            (
                                "genericPointPatchField<Type>::"
                                "genericPointPatchField"
                                "(const pointPatch&, const Field<Type>&, "
                                "const dictionary&)",
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->dimensionedInternalField().name()
                                << " in file "
                                << this->dimensionedInternalField().objectPath()
                                << exit(FatalIOError);
                        }

                        scalarFields_.insert(iter().keyword(), fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<vector> >::typeName
                    )
                    {
                        vectorField* fPtr = new vectorField;
                        fPtr->transfer
                        (
                            dynamicCast<token::Compound<List<vector> > >
                            (
                                fieldToken.transferCompoundToken()
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorIn
                            (
                                "genericPointPatchField<Type>::"
                                "genericPointPatchField"
                                "(const pointPatch&, const Field<Type>&, "
                                "const dictionary&)",
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->dimensionedInternalField().name()
                                << " in file "
                                << this->dimensionedInternalField().objectPath()
                                << exit(FatalIOError);
                        }

                        vectorFields_.insert(iter().keyword(), fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<sphericalTensor> >::typeName
                    )
                    {
                        sphericalTensorField* fPtr = new sphericalTensorField;
                        fPtr->transfer
                        (
                            dynamicCast
                            <
                                token::Compound<List<sphericalTensor> >
                            >
                            (
                                fieldToken.transferCompoundToken()
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorIn
                            (
                                "genericPointPatchField<Type>::"
                                "genericPointPatchField"
                                "(const pointPatch&, const Field<Type>&, "
                                "const dictionary&)",
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->dimensionedInternalField().name()
                                << " in file "
                                << this->dimensionedInternalField().objectPath()
                                << exit(FatalIOError);
                        }

                        sphericalTensorFields_.insert(iter().keyword(), fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<symmTensor> >::typeName
                    )
                    {
                        symmTensorField* fPtr = new symmTensorField;
                        fPtr->transfer
                        (
                            dynamicCast
                            <
                                token::Compound<List<symmTensor> >
                            >
                            (
                                fieldToken.transferCompoundToken()
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorIn
                            (
                                "genericPointPatchField<Type>::"
                                "genericPointPatchField"
                                "(const pointPatch&, const Field<Type>&, "
                                "const dictionary&)",
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->dimensionedInternalField().name()
                                << " in file "
                                << this->dimensionedInternalField().objectPath()
                                << exit(FatalIOError);
                        }

                        symmTensorFields_.insert(iter().keyword(), fPtr);
                    }
                    else if
                    (
                        fieldToken.compoundToken().type()
                     == token::Compound<List<tensor> >::typeName
                    )
                    {
                        tensorField* fPtr = new tensorField;
                        fPtr->transfer
                        (
                            dynamicCast<token::Compound<List<tensor> > >
                            (
                                fieldToken.transferCompoundToken()
                            )
                        );

                        if (fPtr->size() != this->size())
                        {
                            FatalIOErrorIn
                            (
                                "genericPointPatchField<Type>::"
                                "genericPointPatchField"
                                "(const pointPatch&, const Field<Type>&, "
                                "const dictionary&)",
                                dict
                            )   << "\n    size of field " << iter().keyword()
                                << " (" << fPtr->size() << ')'
                                << " is not the same size as the patch ("
                                << this->size() << ')'
                                << "\n    on patch " << this->patch().name()
                                << " of field "
                                << this->dimensionedInternalField().name()
                                << " in file "
                                << this->dimensionedInternalField().objectPath()
                                << exit(FatalIOError);
                        }

                        tensorFields_.insert(iter().keyword(), fPtr);
                    }
                    else
                    {
                        FatalIOErrorIn
                        (
                            "genericPointPatchField<Type>::"
                            "genericPointPatchField"
                            "(const pointPatch&, const Field<Type>&, "
                            "const dictionary&)",
                            dict
                        )   << "\n    compound " << fieldToken.compoundToken()
                            << " not supported"
                            << "\n    on patch " << this->patch().name()
                            << " of field "
                            << this->dimensionedInternalField().name()
                            << " in file "
                            << this->dimensionedInternalField().objectPath()
                            << exit(FatalIOError);
                    }
                }
            }
        }
    }
}


template<class Type>
genericPointPatchField<Type>::genericPointPatchField
(
    const genericPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    calculatedPointPatchField<Type>(ptf, p, iF, mapper),
    actualTypeName_(ptf.actualTypeName_),
    dict_(ptf.dict_)
{
    for
    (
        HashPtrTable<scalarField>::const_iterator iter =
            ptf.scalarFields_.begin();
        iter != ptf.scalarFields_.end();
        ++iter
    )
    {
        scalarFields_.insert(iter.key(), new scalarField(*iter(), mapper));
    }

    for
    (
        HashPtrTable<vectorField>::const_iterator iter =
            ptf.vectorFields_.begin();
        iter != ptf.vectorFields_.end();
        ++iter
    )
    {
        vectorFields_.insert(iter.key(), new vectorField(*iter(), mapper));
    }

    for
    (
        HashPtrTable<sphericalTensorField>::const_iterator iter =
            ptf.sphericalTensorFields_.begin();
        iter != ptf.sphericalTensorFields_.end();
        ++iter
    )
    {
        sphericalTensorFields_.insert
        (
            iter.key(),
            new sphericalTensorField(*iter(), mapper)
        );
    }

    for
    (
        HashPtrTable<symmTensorField>::const_iterator iter =
            ptf.symmTensorFields_.begin();
        iter != ptf.symmTensorFields_.end();
        ++iter
    )
    {
        symmTensorFields_.insert
        (
            iter.key(),
            new symmTensorField(*iter(), mapper)
        );
    }

    for
    (
        HashPtrTable<tensorField>::const_iterator iter =
            ptf.tensorFields_.begin();
        iter != ptf.tensorFields_.end();
        ++iter
    )
    {
        tensorFields_.insert(iter.key(), new tensorField(*iter(), mapper));
    }
}


template<class Type>
genericPointPatchField<Type>::genericPointPatchField
(
    const genericPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    calculatedPointPatchField<Type>(ptf, iF),
    actualTypeName_(ptf.actualTypeName_),
    dict_(ptf.dict_),
    scalarFields_(ptf.scalarFields_),
    vectorFields_(ptf.vectorFields_),
    sphericalTensorFields_(ptf.sphericalTensorFields_),
    symmTensorFields_(ptf.symmTensorFields_),
    tensorFields_(ptf.tensorFields_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void genericPointPatchField<Type>::autoMap
(
    const pointPatchFieldMapper& m
)
{
    for
    (
        HashPtrTable<scalarField>::iterator iter = scalarFields_.begin();
        iter != scalarFields_.end();
        ++iter
    )
    {
        iter()->autoMap(m);
    }

    for
    (
        HashPtrTable<vectorField>::iterator iter = vectorFields_.begin();
        iter != vectorFields_.end();
        ++iter
    )
    {
        iter()->autoMap(m);
    }

    for
    (
        HashPtrTable<sphericalTensorField>::iterator iter =
            sphericalTensorFields_.begin();
        iter != sphericalTensorFields_.end();
        ++iter
    )
    {
        iter()->autoMap(m);
    }

    for
    (
        HashPtrTable<symmTensorField>::iterator iter =
            symmTensorFields_.begin();
        iter != symmTensorFields_.end();
        ++iter
    )
    {
        iter()->autoMap(m);
    }

    for
    (
        HashPtrTable<tensorField>::iterator iter = tensorFields_.begin();
        iter != tensorFields_.end();
        ++iter
    )
    {
        iter()->autoMap(m);
    }
}


template<class Type>
void genericPointPatchField<Type>::rmap
(
    const pointPatchField<Type>& ptf,
    const labelList& addr
)
{
    const genericPointPatchField<Type>& dptf =
        refCast<const genericPointPatchField<Type> >(ptf);

    for
    (
        HashPtrTable<scalarField>::iterator iter = scalarFields_.begin();
        iter != scalarFields_.end();
        ++iter
    )
    {
        HashPtrTable<scalarField>::const_iterator dptfIter =
            dptf.scalarFields_.find(iter.key());

        if (dptfIter != scalarFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    for
    (
        HashPtrTable<vectorField>::iterator iter = vectorFields_.begin();
        iter != vectorFields_.end();
        ++iter
    )
    {
        HashPtrTable<vectorField>::const_iterator dptfIter =
            dptf.vectorFields_.find(iter.key());

        if (dptfIter != vectorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    for
    (
        HashPtrTable<sphericalTensorField>::iterator iter =
            sphericalTensorFields_.begin();
        iter != sphericalTensorFields_.end();
        ++iter
    )
    {
        HashPtrTable<sphericalTensorField>::const_iterator dptfIter =
            dptf.sphericalTensorFields_.find(iter.key());

        if (dptfIter != sphericalTensorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    for
    (
        HashPtrTable<symmTensorField>::iterator iter =
            symmTensorFields_.begin();
        iter != symmTensorFields_.end();
        ++iter
    )
    {
        HashPtrTable<symmTensorField>::const_iterator dptfIter =
            dptf.symmTensorFields_.find(iter.key());

        if (dptfIter != symmTensorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }

    for
    (
        HashPtrTable<tensorField>::iterator iter = tensorFields_.begin();
        iter != tensorFields_.end();
        ++iter
    )
    {
        HashPtrTable<tensorField>::const_iterator dptfIter =
            dptf.tensorFields_.find(iter.key());

        if (dptfIter != tensorFields_.end())
        {
            iter()->rmap(*dptfIter(), addr);
        }
    }
}


template<class Type>
void genericPointPatchField<Type>::write(Ostream& os) const
{
    os.writeKeyword("type") << actualTypeName_ << token::END_STATEMENT << nl;

    for
    (
        dictionary::const_iterator iter = dict_.begin();
        iter != dict_.end();
        ++iter
    )
    {
        if (iter().keyword() != "type")
        {
            if
            (
                iter().isStream()
             && iter().stream().size()
             && iter().stream()[0].isWord()
             && iter().stream()[0].wordToken() == "nonuniform"
            )
            {
                if (scalarFields_.found(iter().keyword()))
                {
                    scalarFields_.find(iter().keyword())()
                        ->writeEntry(iter().keyword(), os);
                }
                else if (vectorFields_.found(iter().keyword()))
                {
                    vectorFields_.find(iter().keyword())()
                        ->writeEntry(iter().keyword(), os);
                }
                else if (sphericalTensorFields_.found(iter().keyword()))
                {
                    sphericalTensorFields_.find(iter().keyword())()
                        ->writeEntry(iter().keyword(), os);
                }
                else if (symmTensorFields_.found(iter().keyword()))
                {
                    symmTensorFields_.find(iter().keyword())()
                        ->writeEntry(iter().keyword(), os);
                }
                else if (tensorFields_.found(iter().keyword()))
                {
                    tensorFields_.find(iter().keyword())()
                        ->writeEntry(iter().keyword(), os);
                }
            }
            else
            {
               iter().write(os);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
