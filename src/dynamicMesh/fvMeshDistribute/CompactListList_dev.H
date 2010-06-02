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

Class
    Foam::CompactListList_dev

Description
    A packed storage unstructured matrix of objects of type \<T\>
    using an offset table for access.

    The offset table is the size of the number of rows+1
    whose elements are the
    accumulated sizes of the rows, i.e.
      - offset[i] gives the index of first element of row i
      - offset[i+1] - offset[i] is the number of elements in row i

    Storage is allocated on free-store during construction.

    As a special case a null-contructed CompactListList_dev has an empty
    offsets_ (instead of size 1).

SourceFiles
    CompactListList_dev.C
    CompactListList_devI.H
    CompactListList_devIO.C

\*---------------------------------------------------------------------------*/

#ifndef CompactListListDev_H
#define CompactListListDev_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T, class Container> class CompactListList_dev;

template<class T, class Container> Istream& operator>>
(
    Istream&,
    CompactListList_dev<T, Container>&
);
template<class T, class Container> Ostream& operator<<
(
    Ostream&,
    const CompactListList_dev<T, Container>&
);


/*---------------------------------------------------------------------------*\
                       Class CompactListList_dev Declaration
\*---------------------------------------------------------------------------*/

template<class T, class Container = List<T> >
class CompactListList_dev
{
    // Private data

        label size_;

        //- Offset table
        List<label> offsets_;

        //- Packed matrix of data
        List<T> m_;


public:

    // Static Member Functions

        //- Return a null CompactListList_dev
        inline static const CompactListList_dev<T, Container>& null();

    // Constructors

        //- Null constructor.
        inline CompactListList_dev();

        //- Construct by converting given List<List<T> >
        explicit CompactListList_dev(const List<Container>&);

        //- Construct given size of offset table (number of rows)
        //  and number of data.
        inline CompactListList_dev(const label nRows, const label nData);

        //- Construct given size of offset table (number of rows),
        //  the number of data and a value for all elements.
        inline CompactListList_dev(const label nRows, const label nData, const T&);

        //- Construct given list of row-sizes.
        explicit CompactListList_dev(const UList<label>& rowSizes);

        //- Construct given list of row-sizes
        CompactListList_dev(const UList<label>& rowSizes, const T&);

        //- Construct by transferring the parameter contents
        explicit CompactListList_dev(const Xfer<CompactListList_dev<T, Container> >&);

        //- Construct as copy or re-use as specified.
        CompactListList_dev(CompactListList_dev<T, Container>&, bool reUse);

        //- Construct from Istream.
        CompactListList_dev(Istream&);

        //- Clone
        inline autoPtr<CompactListList_dev<T, Container> > clone() const;


    // Member Functions

        // Access

            //- Return the primary size, i.e. the number of rows
            inline label size() const;

            //- Return true if the number of rows is zero
            inline bool empty() const;

            //- Return the offset table (= size()+1)
            inline const List<label>& offsets() const;

            //- Return non-const access to the offset table
            inline List<label>& offsets();

            //- Return the packed matrix of data
            inline const List<T>& m() const;

            //- Return non-const access to the packed matrix of data
            inline List<T>& m();


        // Edit

            //- Reset size of CompactListList_dev.
            //  This form only allows contraction of the CompactListList_dev.
            void setSize(const label nRows);

            //- Reset size of CompactListList_dev.
            void setSize(const label nRows, const label nData);

            //- Reset sizes of CompactListList_dev and value for new elements.
            void setSize(const label nRows, const label nData, const T&);

            //- Reset size of CompactListList_dev.
            void setSize(const UList<label>& rowSizes);

            //- Reset size of CompactListList_dev.
            //  This form only allows contraction of the CompactListList_dev.
            inline void resize(const label nRows);

            //- Reset size of CompactListList_dev.
            inline void resize(const label nRows, const label nData);

            //- Reset sizes of CompactListList_dev and value for new elements.
            inline void resize(const label nRows, const label nData, const T&);

            //- Reset size of CompactListList_dev.
            inline void resize(const UList<label>& rowSizes);

            //- Clear the CompactListList_dev, i.e. set sizes to zero.
            void clear();

            //- Return sizes (to be used e.g. for construction)
            labelList sizes() const;

            //- Transfer the contents of the argument CompactListList_dev
            //  into this CompactListList_dev and annul the argument list.
            void transfer(CompactListList_dev<T, Container>&);

            //- Transfer the contents to the Xfer container
            inline Xfer<CompactListList_dev<T, Container> > xfer();

        // Other

            //- Return index into m
            inline label index(const label row, const label col) const;

            //- Get row for index into m.
            inline label whichRow(const label index) const;

            //- Get column index (j) given above row
            inline label whichColumn(const label row, const label index) const;


    // Member operators

        //- Return subscript-checked row as UList.
        inline UList<T> operator[](const label i);

        //- Return const subscript-checked row as UList.
        inline const UList<T> operator[](const label i) const;

        //- Return subscript-checked element.
        inline T& operator()(const label i, const label j);

        //- Return const subscript-checked element.
        inline const T& operator()(const label i, const label j) const;

        //- Return as List<Container>
        List<Container> operator()() const;

        //- Assignment of all entries to the given value
        inline void operator=(const T&);


    // Istream operator

        //- Read CompactListList_dev from Istream, discarding contents
        //  of existing CompactListList_dev.
        friend Istream& operator>> <T, Container>
        (
            Istream&,
            CompactListList_dev<T, Container>&
        );

        // Write CompactListList_dev to Ostream.
        friend Ostream& operator<< <T, Container>
        (
            Ostream&,
            const CompactListList_dev<T, Container>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "CompactListList_devI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CompactListList_dev.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
