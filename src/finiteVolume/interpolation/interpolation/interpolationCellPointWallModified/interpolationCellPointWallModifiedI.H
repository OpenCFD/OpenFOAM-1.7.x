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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::interpolationCellPointWallModified<Type>::interpolate
(
    const cellPointWeightWallModified& cpw
) const
{
    const FixedList<scalar, 4>& weights = cpw.weights();
    const FixedList<label, 3>& faceVertices = cpw.faceVertices();

    Type t = this->psip_[faceVertices[0]]*weights[0];
    t += this->psip_[faceVertices[1]]*weights[1];
    t += this->psip_[faceVertices[2]]*weights[2];
    t += this->psi_[cpw.cell()]*weights[3];

    return t;
}


template<class Type>
inline Type Foam::interpolationCellPointWallModified<Type>::interpolate
(
    const vector& position,
    const label celli,
    const label facei
) const
{
    return
        interpolate
        (
            cellPointWeightWallModified
            (
                this->pMesh_,
                position,
                celli,
                facei
            )
        );
}


// ************************************************************************* //
