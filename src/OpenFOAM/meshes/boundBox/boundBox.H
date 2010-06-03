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
    Foam::boundBox

Description
    A bounding box defined in terms of the points at its extremities.

\*---------------------------------------------------------------------------*/

#ifndef boundBox_H
#define boundBox_H

#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class boundBox;
template<class T> class tmp;

Ostream& operator<<(Ostream& os, const boundBox& b);


/*---------------------------------------------------------------------------*\
                          Class boundBox Declaration
\*---------------------------------------------------------------------------*/

class boundBox
{
    // Private data

        //- Minimum and maximum describing the bounding box
        point min_, max_;

    // Private Member Functions

        //- Calculate the bounding box from the given pointField.
        //  Does parallel communication (doReduce = true)
        void calculate(const pointField&, const bool doReduce = true);

public:

    // Static data members

        //- The great value used for greatBox and invertedBox
        static const scalar great;

        //- A very large boundBox: min/max == -/+ VGREAT
        static const boundBox greatBox;

        //- A very large inverted boundBox: min/max == +/- VGREAT
        static const boundBox invertedBox;


    // Constructors

        //- Construct null, setting points to zero
        boundBox()
        :
            min_(point::zero),
            max_(point::zero)
        {}

        //- Construct from components
        boundBox(const point& min, const point& max)
        :
            min_(min),
            max_(max)
        {}

        //- Construct as the bounding box of the given pointField.
        //  Does parallel communication (doReduce = true)
        boundBox(const pointField&, const bool doReduce = true);

        //- Construct as the bounding box of the given temporary pointField.
        //  Does parallel communication (doReduce = true)
        boundBox(const tmp<pointField>&, const bool doReduce = true);

        //- Construct from Istream
        boundBox(Istream&);


    // Member functions

        // Access

            //- Minimum describing the bounding box
            const point& min() const
            {
                return min_;
            }

            //- Maximum describing the bounding box
            const point& max() const
            {
                return max_;
            }

            //- Minimum describing the bounding box, non-const access
            point& min()
            {
                return min_;
            }

            //- Maximum describing the bounding box, non-const access
            point& max()
            {
                return max_;
            }

            //- The midpoint of the bounding box
            point midpoint() const
            {
                return 0.5 * (max_ + min_);
            }

            //- The bounding box span (from minimum to maximum)
            vector span() const
            {
                return (max_ - min_);
            }

            //- The magnitude of the bounding box span
            scalar mag() const
            {
                return ::Foam::mag(max_ - min_);
            }

            //- Smallest length/height/width dimension
            scalar minDim() const
            {
                return cmptMin(span());
            }

            //- Largest length/height/width dimension
            scalar maxDim() const
            {
                return cmptMax(span());
            }

            //- Average length/height/width dimension
            scalar avgDim() const
            {
                return cmptAv(span());
            }


        // Query

            //- Overlaps/touches boundingBox?
            bool overlaps(const boundBox& bb) const
            {
                return
                (
                    bb.max_.x() >= min_.x() && bb.min_.x() <= max_.x()
                 && bb.max_.y() >= min_.y() && bb.min_.y() <= max_.y()
                 && bb.max_.z() >= min_.z() && bb.min_.z() <= max_.z()
                );
            }

            //- Contains point? (inside or on edge)
            bool contains(const point& pt) const
            {
                return
                (
                    pt.x() >= min_.x() && pt.x() <= max_.x()
                 && pt.y() >= min_.y() && pt.y() <= max_.y()
                 && pt.z() >= min_.z() && pt.z() <= max_.z()
                );
            }

            //- Contains point? (inside only)
            bool containsInside(const point& pt) const
            {
                return
                (
                    pt.x() > min_.x() && pt.x() < max_.x()
                 && pt.y() > min_.y() && pt.y() < max_.y()
                 && pt.z() > min_.z() && pt.z() < max_.z()
                );
            }


    // Friend Operators

        friend bool operator==(const boundBox& a, const boundBox& b)
        {
            return (a.min_ == b.min_) && (a.max_ == b.max_);
        }

        friend bool operator!=(const boundBox& a, const boundBox& b)
        {
            return !(a == b);
        }


    // IOstream operator

        friend Istream& operator>>(Istream& is, boundBox&);
        friend Ostream& operator<<(Ostream& os, const boundBox&);
};


//- Data associated with boundBox type are contiguous
template<>
inline bool contiguous<boundBox>() {return contiguous<point>();}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// #include "boundBoxI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
