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

#include "interpolatePointToCell.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::interpolatePointToCell
(
    const GeometricField<Type, pointPatchField, pointMesh>& ptf,
    const label celli
)
{
    const primitiveMesh& mesh = ptf.mesh()();

    const cell& cFaces = mesh.cells()[celli];

    labelHashSet pointHad(10*cFaces.size());

    Type sum(pTraits<Type>::zero);

    forAll(cFaces, i)
    {
        const face& f = mesh.faces()[cFaces[i]];

        forAll(f, fp)
        {
            label v = f[fp];

            if (pointHad.insert(v))
            {
                sum += ptf[v];
            }
        }
    }

    return sum/pointHad.size();
}


// ************************************************************************* //
