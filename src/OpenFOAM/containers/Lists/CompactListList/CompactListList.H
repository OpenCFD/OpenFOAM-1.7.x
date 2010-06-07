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
    Foam::CompactListList

Description
    A packed storage unstructured matrix of objects of type \<T\>
    using an offset table for access.

    The offset table is the size of the number of rows whose elements are the
    accumulated sizes of the rows, i.e.
      - offset[i] gives the index of first element of row i + 1
      - offset[i] - offset[i-1] is the number of elements in row i

    and for i = 0, offset[i-1] = 0.

    Storage is allocated on free-store during construction.

SourceFiles
    CompactListList.C
    CompactListListI.H
    CompactListListIO.C

\*---------------------------------------------------------------------------*/

#ifndef CompactListList_H
#define CompactListList_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T> class CompactListList;

template<class T> Istream& operator>>(Istream&, CompactListList<T>&);
template<class T> Ostream& operator<<(Ostream&, const CompactListList<T>&);


/*---------------------------------------------------------------------------*\
                       Class CompactListList Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class CompactListList
{
    // Private data

        //- Offset table
        List<label> offsets_;

        //- Packed matrix of data
        List<T> m_;


public:

    // Static Member Functions

        //- Return a null CompactListList
        inline static const CompactListList<T>& null();

    // Constructors

        //- Null constructor.
        inline CompactListList();

        //- Construct by converting given List<List<T> >
        CompactListList(const List<List<T> >&);

        //- Construct given size of offset table (number of rows)
        //  and number of data.
        inline CompactListList(const label nRows, const label nData);

        //- Construct given size of offset table (number of rows),
        //  the number of data and a value for all elements.
        inline CompactListList(const label nRows, const label nData, const T&);

        //- Construct given list of row-sizes.
        CompactListList(const UList<label>& rowSizes);

        //- Construct given list of row-sizes
        CompactListList(const UList<label>& rowSizes, const T&);

        //- Construct by transferring the parameter contents
        CompactListList(const Xfer<CompactListList<T> >&);

        //- Construct as copy or re-use as specified.
        CompactListList(CompactListList<T>&, bool reUse);

        //- Construct from Istream.
        CompactListList(Istream&);

        //- Clone
        inline autoPtr<CompactListList<T> > clone() const;


    // Member Functions

        // Access

            //- Return the primary size, i.e. the number of rows
            inline label size() const;

            //- Return true if the number of rows is zero
            inline bool empty() const;

            //- Return the offset table
            inline const List<label>& offsets() const;

            //- Return non-const access to the offset table
            inline List<label>& offsets();

            //- Return the packed matrix of data
            inline const List<T>& m() const;

            //- Return non-const access to the packed matrix of data
            inline List<T>& m();


        // Edit

            //- Reset size of CompactListList.
            //  This form only allows contraction of the CompactListList.
            void setSize(const label nRows);

            //- Reset size of CompactListList.
            void setSize(const label nRows, const label nData);

            //- Reset sizes of CompactListList and value for new elements.
            void setSize(const label nRows, const label nData, const T&);

            //- Reset size of CompactListList.
            void setSize(const UList<label>& rowSizes);

            //- Reset size of CompactListList.
            //  This form only allows contraction of the CompactListList.
            inline void resize(const label nRows);

            //- Reset size of CompactListList.
            inline void resize(const label nRows, const label nData);

            //- Reset sizes of CompactListList and value for new elements.
            inline void resize(const label nRows, const label nData, const T&);

            //- Reset size of CompactListList.
            inline void resize(const UList<label>& rowSizes);

            //- Clear the CompactListList, i.e. set sizes to zero.
            void clear();

            //- Return sizes (to be used e.g. for construction)
            labelList sizes() const;

            //- Transfer the contents of the argument CompactListList
            //  into this CompactListList and annull the argument list.
            void transfer(CompactListList<T>&);

            //- Transfer the contents to the Xfer container
            inline Xfer<CompactListList<T> > xfer();

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

        //- Return as List<List<T> >
        List<List<T> > operator()() const;

        //- Assignment of all entries to the given value
        inline void operator=(const T&);


    // Istream operator

        //- Read CompactListList from Istream, discarding contents
        //  of existing CompactListList.
        friend Istream& operator>> <T>(Istream&, CompactListList<T>&);

        // Write CompactListList to Ostream.
        friend Ostream& operator<< <T>(Ostream&, const CompactListList<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "CompactListListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CompactListList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
