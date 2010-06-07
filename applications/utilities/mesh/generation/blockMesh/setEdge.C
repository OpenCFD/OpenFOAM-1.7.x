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

Description
    from the list of curved edges creates a list
    of edges that are not curved. It is assumed
    that all other edges are straight lines

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "blockDescriptor.H"
#include "lineEdge.H"
#include "lineDivide.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar calcGexp(const scalar expRatio, const label dim)
{
    if (dim == 1)
    {
        return 0.0;
    }
    else
    {
        return pow(expRatio, 1.0/(dim - 1));
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void blockDescriptor::setEdge(label edgeI, label start, label end, label dim)
{
    // for all edges check the list of curved edges. If the edge is curved,
    // add it to the list. If the edge is not found, create is as a line

    bool found = false;

    // set reference to the list of labels defining the block
    const labelList& blockLabels = blockShape_;

    // set reference to global list of points
    const pointField blockPoints = blockShape_.points(blockMeshPoints_);

    // x1
    found = false;

    forAll (curvedEdges_, nCEI)
    {
        if (curvedEdges_[nCEI].compare(blockLabels[start], blockLabels[end]))
        {
            found = true;

            // check the orientation:
            // if the starting point of the curve is the same as the starting
            // point of the edge, do the parametrisation and pick up the points
            if (blockLabels[start] == curvedEdges_[nCEI].start())
            {
                // calculate the geometric expension factor out of the
                // expansion ratio
                scalar gExp = calcGexp(expand_[edgeI], dim);

                // divide the line
                lineDivide divEdge(curvedEdges_[nCEI], dim, gExp);

                edgePoints_[edgeI] = divEdge.points();
                edgeWeights_[edgeI] = divEdge.lambdaDivisions();
            }
            else
            {
                // the curve has got the opposite orientation
                scalar gExp = calcGexp(expand_[edgeI], dim);

                // divide the line
                lineDivide divEdge(curvedEdges_[nCEI], dim, 1.0/(gExp+SMALL));

                pointField p = divEdge.points();
                scalarList d = divEdge.lambdaDivisions();

                edgePoints_[edgeI].setSize(p.size());
                edgeWeights_[edgeI].setSize(d.size());

                label pMax = p.size() - 1;
                forAll (p, pI)
                {
                    edgePoints_[edgeI][pI] = p[pMax - pI];
                    edgeWeights_[edgeI][pI] = 1.0 - d[pMax - pI];
                }
            }

            break;
        }
    }

    if (!found)
    {
        // edge is a straight line
        scalar gExp = calcGexp(expand_[edgeI], dim);

        // divide the line
        lineDivide divEdge
        (
            lineEdge(blockPoints, start, end),
            dim,
            gExp
        );

        edgePoints_[edgeI] = divEdge.points();
        edgeWeights_[edgeI] = divEdge.lambdaDivisions();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
