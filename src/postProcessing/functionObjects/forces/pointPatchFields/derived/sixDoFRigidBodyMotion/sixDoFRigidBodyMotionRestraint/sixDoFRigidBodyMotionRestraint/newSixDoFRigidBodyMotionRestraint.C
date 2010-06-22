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

#include "sixDoFRigidBodyMotionRestraint.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::sixDoFRigidBodyMotionRestraint>
Foam::sixDoFRigidBodyMotionRestraint::New(const dictionary& sDoFRBMRDict)
{
    word sixDoFRigidBodyMotionRestraintTypeName =
        sDoFRBMRDict.lookup("sixDoFRigidBodyMotionRestraint");

    // Info<< "Selecting sixDoFRigidBodyMotionRestraint function "
    //     << sixDoFRigidBodyMotionRestraintTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
    dictionaryConstructorTablePtr_->find
    (
        sixDoFRigidBodyMotionRestraintTypeName
    );

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "sixDoFRigidBodyMotionRestraint::New"
            "("
                "const dictionary& sDoFRBMRDict"
            ")"
        )   << "Unknown sixDoFRigidBodyMotionRestraint type "
            << sixDoFRigidBodyMotionRestraintTypeName << endl << endl
            << "Valid  sixDoFRigidBodyMotionRestraints are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<sixDoFRigidBodyMotionRestraint>(cstrIter()(sDoFRBMRDict));
}


// ************************************************************************* //
