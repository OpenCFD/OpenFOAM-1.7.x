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

#include "Kmesh.H"
#include "polyMesh.H"
#include "volFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   //! @cond fileScope
   inline label rep
   (
       const label i,
       const label j,
       const label k,
       const labelList& nn
   )
   {
       return (k + j*nn[2] + i*nn[1]*nn[2]);
   }
   //! @endcond fileScope

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// from fvMesh
Foam::Kmesh::Kmesh(const fvMesh& mesh)
:
    vectorField(mesh.V().size()),
    NN(vector::dim)
{
    const scalar pi = mathematicalConstant::pi;
    const scalar twoPi = 2.0*pi;

    boundBox box = mesh.bounds();
    L = box.span();

    vector cornerCellCentre = ::Foam::max(mesh.C().internalField());
    vector cellL = 2 * (box.max() - cornerCellCentre);

    vector rdeltaByL;
    label nTot = 1;

    label i;
    forAll(NN, i)
    {
        NN[i] = label(L[i]/cellL[i] + 0.5);
        nTot *= NN[i];

        if (NN[i] > 1)
        {
            L[i] -= cellL[i];
        }

        rdeltaByL[i] = NN[i]/(L[i]*L[i]);
    }

    if (nTot != mesh.nCells())
    {
        FatalErrorIn("Kmesh::Kmesh(const fvMesh& mesh)")
            << "calculated number of cells is incorrect"
            << abort(FatalError);
    }

    for (i=0; i<NN[0]; i++)
    {
        scalar k1 = (i - NN[0]/2)*twoPi/L[0];

        for (label j=0; j<NN[1]; j++)
        {
            scalar k2 = (j - NN[1]/2)*twoPi/L[1];

            for (label k=0; k<NN[2]; k++)
            {
                scalar k3 = (k - NN[2]/2)*twoPi/L[2];

                (*this)[rep(i, j, k, NN)] = vector(k1, k2, k3);
            }
        }
    }

    Kmax = mag((*this)[rep(NN[0]-1, NN[1]-1, NN[2]-1, NN)]);
}


// ************************************************************************* //
