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

#include "PackedList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<unsigned nBits>
Foam::PackedList<nBits>::PackedList(const label size, const unsigned int val)
:
    StorageList(packedLength(size), 0u),
    size_(size)
{
    operator=(val);
}


template<unsigned nBits>
Foam::PackedList<nBits>::PackedList(const UList<label>& lst)
:
    StorageList(packedLength(lst.size()), 0u),
    size_(lst.size())
{
    forAll(lst, i)
    {
        set(i, lst[i]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


#if (UINT_MAX == 0xFFFFFFFF)
// 32-bit counting, Hamming weight method
#   define COUNT_PACKEDBITS(sum, x)                                           \
{                                                                             \
    x -= (x >> 1) & 0x55555555;                                               \
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);                           \
    sum += (((x + (x >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;                \
}
#elif (UINT_MAX == 0xFFFFFFFFFFFFFFFF)
// 64-bit counting, Hamming weight method
#   define COUNT_PACKEDBITS(sum, x)                                           \
{                                                                             \
    x -= (x >> 1) & 0x5555555555555555;                                       \
    x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);           \
    sum += (((x + (x >> 4)) & 0x0F0F0F0F0F0F0F0F) * 0x0101010101010101) >> 56;\
}
#else
// Arbitrary number of bits, Brian Kernighan's method
#   define COUNT_PACKEDBITS(sum, x)    for (; x; ++sum) { x &= x - 1; }
#endif


template<unsigned nBits>
unsigned int Foam::PackedList<nBits>::count() const
{
    register unsigned int c = 0;

    if (size_)
    {
        // mask value for complete segments
        unsigned int mask = maskLower(packing());

        const unsigned int endSeg = size_ / packing();
        const unsigned int endOff = size_ % packing();

        // count bits in complete segments
        for (unsigned i = 0; i < endSeg; ++i)
        {
            register unsigned int bits = StorageList::operator[](i) & mask;
            COUNT_PACKEDBITS(c, bits);
        }

        // count bits in partial segment
        if (endOff)
        {
            mask = maskLower(endOff);

            register unsigned int bits = StorageList::operator[](endSeg) & mask;
            COUNT_PACKEDBITS(c, bits);
        }
    }

    return c;
}


template<unsigned nBits>
bool Foam::PackedList<nBits>::trim()
{
    if (!size_)
    {
        return false;
    }

    // mask value for complete segments
    unsigned int mask = maskLower(packing());

    label currElem = packedLength(size_) - 1;
    unsigned int endOff = size_ % packing();

    // clear trailing bits on final segment
    if (endOff)
    {
        StorageList::operator[](currElem) &= maskLower(endOff);
    }

    // test entire segment
    while (currElem > 0 && !(StorageList::operator[](currElem) &= mask))
    {
        currElem--;
    }

    // test segment
    label newsize = (currElem + 1) * packing();

    // mask for the final segment
    mask = max_value() << (nBits * (packing() - 1));

    for (endOff = packing(); endOff >= 1; --endOff, --newsize)
    {
        if (StorageList::operator[](currElem) & mask)
        {
            break;
        }

        mask >>= nBits;
    }

    if (size_ == newsize)
    {
        return false;
    }

    size_ = newsize;
    return false;
}


template<unsigned nBits>
void Foam::PackedList<nBits>::flip()
{
    label packLen = packedLength(size_);

    for (label i=0; i < packLen; i++)
    {
        StorageList::operator[](i) = ~StorageList::operator[](i);
    }
}


template<unsigned nBits>
Foam::labelList Foam::PackedList<nBits>::values() const
{
    labelList elems(size_);

    forAll(*this, i)
    {
        elems[i] = get(i);
    }
    return elems;
}


template<unsigned nBits>
Foam::Ostream& Foam::PackedList<nBits>::iteratorBase::print(Ostream& os) const
{
    os  << "iterator<"  << label(nBits) << "> ["
        << this->index_ << "]"
        << " segment:"  << label(this->index_ / packing())
        << " offset:"   << label(this->index_ % packing())
        << " value:"    << this->get()
        << nl;

    return os;
}


template<unsigned nBits>
Foam::Ostream& Foam::PackedList<nBits>::print(Ostream& os) const
{
    const label packLen = packedLength(size_);

    os  << "PackedList<" << nBits << ">"
        << " max_value:" << max_value()
        << " packing:"   << packing() << nl
        << " count: "     << count() << nl
        << " size/capacity: " << size_ << "/" << capacity() << nl
        << " storage/capacity: " << packLen << "/" << StorageList::size()
        << "\n(\n";

    // mask value for complete segments
    unsigned int mask = maskLower(packing());

    for (label i=0; i < packLen; i++)
    {
        const StorageType& rawBits = StorageList::operator[](i);

        // the final segment may not be full, modify mask accordingly
        if (i+1 == packLen)
        {
            unsigned int endOff = size_ % packing();

            if (endOff)
            {
                mask = maskLower(endOff);
            }
            else
            {
                continue;
            }
        }

        for (unsigned int testBit = (1u << max_bits()); testBit; testBit >>= 1)
        {
            if (mask & testBit)
            {
                if (rawBits & testBit)
                {
                    os << '1';
                }
                else
                {
                    os << '-';
                }
            }
            else
            {
                os << 'x';
            }
        }
        os << '\n';
    }
    os << ")\n";

    return os;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<unsigned nBits>
void Foam::PackedList<nBits>::operator=(const PackedList<nBits>& lst)
{
    StorageList::operator=(lst);
    size_ = lst.size();
}


template<unsigned nBits>
void Foam::PackedList<nBits>::operator=(const UList<label>& lst)
{
    setCapacity(lst.size());
    size_ = lst.size();

    forAll(lst, i)
    {
        set(i, lst[i]);
    }
}


// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

//template<unsigned nBits>
//Foam::Ostream& ::Foam::operator<<(Ostream& os, const PackedList<nBits>& lst)
//{
//    os << lst();
//    return os;
//}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
