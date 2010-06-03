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

\*---------------------------------------------------------------------------*/

#include "meshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void meshToMesh::calculateInverseDistanceWeights() const
{
    if (debug)
    {
        Info<< "meshToMesh::calculateInverseDistanceWeights() : "
            << "calculating inverse distance weighting factors" << endl;
    }

    if (inverseDistanceWeightsPtr_)
    {
        FatalErrorIn("meshToMesh::calculateInverseDistanceWeights()")
            << "weighting factors already calculated"
            << exit(FatalError);
    }

    inverseDistanceWeightsPtr_ = new scalarListList(toMesh_.nCells());
    scalarListList& invDistCoeffs = *inverseDistanceWeightsPtr_;

    // get reference to source mesh data
    const labelListList& cc = fromMesh_.cellCells();
    const vectorField& centreFrom = fromMesh_.C().internalField();
    const vectorField& centreTo = toMesh_.C().internalField();

    forAll (cellAddressing_, celli)
    {
        if (cellAddressing_[celli] != -1)
        {
            const vector& target = centreTo[celli];
            scalar m = mag(target - centreFrom[cellAddressing_[celli]]);

            const labelList& neighbours = cc[cellAddressing_[celli]];

            // if the nearest cell is a boundary cell or there is a direct hit,
            // pick up the value
            if
            (
                m < directHitTol                            // Direct hit
             || neighbours.empty()
            )
            {
                invDistCoeffs[celli].setSize(1);
                invDistCoeffs[celli][0] = 1.0;
            }
            else
            {
                invDistCoeffs[celli].setSize(neighbours.size() + 1);

                // The first coefficient corresponds to the centre cell.
                // The rest is ordered in the same way as the cellCells list.
                scalar invDist = 1.0/m;
                invDistCoeffs[celli][0] = invDist;
                scalar sumInvDist = invDist;

                // now add the neighbours
                forAll (neighbours, ni)
                {
                    invDist = 1.0/mag(target - centreFrom[neighbours[ni]]);
                    invDistCoeffs[celli][ni + 1] = invDist;
                    sumInvDist += invDist;
                }
                
                // divide by the total inverse-distance
                forAll (invDistCoeffs[celli], i)
                {
                    invDistCoeffs[celli][i] /= sumInvDist;
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const scalarListList& meshToMesh::inverseDistanceWeights() const
{
    if (!inverseDistanceWeightsPtr_)
    {
        calculateInverseDistanceWeights();
    }
    
    return *inverseDistanceWeightsPtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
