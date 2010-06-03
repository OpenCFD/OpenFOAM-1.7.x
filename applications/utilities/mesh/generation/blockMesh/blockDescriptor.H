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
    Foam::blockDescriptor

Description
    block descriptor. Takes the description of the block and the list
    of curved edges and creates a list of points on edges together
    with the weighting factors

SourceFiles
    blockDescriptor.C

\*---------------------------------------------------------------------------*/

#ifndef blockDescriptor_H
#define blockDescriptor_H

#include "scalar.H"
#include "label.H"
#include "point.H"
#include "cellShape.H"
#include "scalarList.H"
#include "curvedEdgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                           Class blockDescriptor Declaration
\*---------------------------------------------------------------------------*/

class blockDescriptor
{
    // Private data

        //- Block mesh points
        const pointField& blockMeshPoints_;

        //- block shape
        cellShape blockShape_;

        // reference to a list of curved edges
        const curvedEdgeList& curvedEdges_;

        // block edge points
        List<List<point> > edgePoints_;

        //- block edge weighting factors
        scalarListList edgeWeights_;

        //- number of point in each direction
        Vector<label> n_;

        //- expansion ratios in all directions
        scalarList expand_;

        //- name of the zone (empty string if none)
        word zoneName_;

    // Private member functions

        void makeBlockEdges();

        void setEdge(label edge, label start, label end, label dim);


public:

    // Constructors

        //- Construct from components. Optional cellSet/zone name.
        blockDescriptor
        (
            const cellShape&,
            const pointField& blockMeshPoints,
            const curvedEdgeList& edges,
            const Vector<label>& n,
            const scalarList& expand,
            const word& zoneName = ""
        );

        //- Construct from Istream
        blockDescriptor
        (
            const pointField& blockMeshPoints,
            const curvedEdgeList& edges,
            Istream& is
        );

        //- Clone
        autoPtr<blockDescriptor> clone() const
        {
            notImplemented("blockDescriptor::clone()");
            return autoPtr<blockDescriptor>(NULL);
        }


    // Member Functions

        // Access

        const pointField& points() const;

        const cellShape& blockShape() const;

        const List<List<point> >& blockEdgePoints() const;

        const scalarListList& blockEdgeWeights() const;

        const Vector<label>& n() const;

        const word& zoneName() const;

    // Member Operators

        void operator=(const blockDescriptor&);


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const blockDescriptor&);
};


inline Istream& operator>>(Istream& is, blockDescriptor*)
{
    notImplemented("Istream& operator>>(Istream& is, blockDescriptor*)");
    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
