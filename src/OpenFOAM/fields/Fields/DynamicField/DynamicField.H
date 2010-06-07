/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::DynamicField

Description
    Dynamically sized Field. WIP.

SourceFiles
    DynamicField.C

\*---------------------------------------------------------------------------*/

#ifndef DynamicField_H
#define DynamicField_H

#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type>
class DynamicField;

template<class Type>
Ostream& operator<<(Ostream&, const DynamicField<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const tmp<DynamicField<Type> >&);

template<class Type>
Istream& operator>>(Istream&, DynamicField<Type>&);


/*---------------------------------------------------------------------------*\
                           Class DynamicField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DynamicField
:
    public Field<Type>
{

    // Private data

        //- The capacity (allocated size) of the underlying field.
        label capacity_;


        //- Construct given size and initial value
        DynamicField(const label, const Type&);

        //- Construct as copy of tmp<DynamicField>
#       ifdef ConstructFromTmp
        DynamicField(const tmp<DynamicField<Type> >&);
#       endif

        //- Construct from a dictionary entry
        DynamicField(const word&, const dictionary&, const label);

public:

    // Static data members

        static const char* const typeName;


    // Static Member Functions

        //- Return a null field
        inline static const DynamicField<Type>& null()
        {
            return *reinterpret_cast< DynamicField<Type>* >(0);
        }


    // Constructors

        //- Construct null
        //  Used for temporary fields which are initialised after construction
        DynamicField();

        //- Construct given size
        //  Used for temporary fields which are initialised after construction
        explicit inline DynamicField(const label);

        //- Construct as copy of a UList\<Type\>
        explicit inline DynamicField(const UList<Type>&);

        //- Construct by transferring the List contents
        explicit inline DynamicField(const Xfer<List<Type> >&);

        //- Construct by 1 to 1 mapping from the given field
        inline DynamicField
        (
            const UList<Type>& mapF,
            const labelList& mapAddressing
        );

        //- Construct by interpolative mapping from the given field
        inline DynamicField
        (
            const UList<Type>& mapF,
            const labelListList& mapAddressing,
            const scalarListList& weights
        );

        //- Construct by mapping from the given field
        inline DynamicField
        (
            const UList<Type>& mapF,
            const FieldMapper& map
        );

        //- Construct as copy
        inline DynamicField(const DynamicField<Type>&);

        //- Construct as copy or re-use as specified.
        inline DynamicField(DynamicField<Type>&, bool reUse);

        //- Construct by transferring the Field contents
        inline DynamicField(const Xfer<DynamicField<Type> >&);

        //- Construct from Istream
        inline DynamicField(Istream&);

        //- Clone
        tmp<DynamicField<Type> > clone() const;


    // Member Functions

        //- Size of the underlying storage.
        inline label capacity() const;

        //- Append an element at the end of the list
        inline void append(const Type&);

        //- Alter the addressed list size.
        //  New space will be allocated if required.
        //  Use this to resize the list prior to using the operator[] for
        //  setting values (as per List usage).
        void setSize(const label nElem);

    // Member operators

        inline void operator=(const DynamicField<Type>&);
        inline void operator=(const UList<Type>&);
        inline void operator=(const tmp<DynamicField<Type> >&);

        //- Return element of Field.
        inline Type& operator[](const label i);

        //- Return element of constant Field.
        inline const Type& operator[](const label) const;

    // IOstream operators

        friend Ostream& operator<< <Type>
        (Ostream&, const DynamicField<Type>&);

        friend Ostream& operator<< <Type>
        (Ostream&, const tmp<DynamicField<Type> >&);

        friend Istream& operator>> <Type>
        (Istream&, DynamicField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DynamicFieldI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DynamicField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
