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

#include "fvMotionSolver.H"
#include "fixedValuePointPatchFields.H"
#include "cellMotionFvPatchFields.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::fvMotionSolver::cellMotionBoundaryTypes
(
    const typename GeometricField<Type, pointPatchField, pointMesh>::
    GeometricBoundaryField& pmUbf
) const
{
    wordList cmUbf = pmUbf.types();

    // Remove global patches from the end of the list
    cmUbf.setSize(fvMesh_.boundary().size());

    forAll(cmUbf, patchi)
    {
        if (isA<fixedValuePointPatchField<Type> >(pmUbf[patchi]))
        {
            cmUbf[patchi] = cellMotionFvPatchField<Type>::typeName;
        }

        if (debug)
        {
            Pout<< "Patch:" << fvMesh_.boundary()[patchi].patch().name()
                << " pointType:" << pmUbf.types()[patchi]
                << " cellType:" << cmUbf[patchi] << endl;
        }
    }

    return cmUbf;
}


// ************************************************************************* //
