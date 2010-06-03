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

#include "error.H"

#include "ORourkeCollisionModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ORourkeCollisionModel, 0);

addToRunTimeSelectionTable
(
    collisionModel,
    ORourkeCollisionModel,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
ORourkeCollisionModel::ORourkeCollisionModel
(
    const dictionary& dict,
    spray& sm,
    Random& rndGen
)
:
    collisionModel(dict, sm, rndGen),
    vols_(sm.mesh().V()),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    coalescence_(coeffsDict_.lookup("coalescence"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

ORourkeCollisionModel::~ORourkeCollisionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ORourkeCollisionModel::collideParcels(const scalar dt) const
{

    if (spray_.size() < 2)
    {
        return;
    }

    spray::iterator secondParcel = spray_.begin();
    ++secondParcel;
    spray::iterator p1 = secondParcel;

    while (p1 != spray_.end())
    {
        label cell1 = p1().cell();

        spray::iterator p2 = spray_.begin();

        while (p2 != p1)
        {
            label cell2 = p2().cell();

            // No collision if parcels are not in the same cell
            if (cell1 == cell2)
            {
#               include "sameCell.H"
            } // if - parcels in same cell

            // remove coalesced droplet
            if (p2().m() < VSMALL) 
            {
                spray::iterator tmpElmnt = p2;
                ++tmpElmnt;
                spray_.deleteParticle(p2());
                p2 = tmpElmnt;
            }
            else 
            {
                ++p2;
            }

        } // inner loop

        // remove coalesced droplet
        if (p1().m() < VSMALL) 
        {
            spray::iterator tmpElmnt = p1;
            ++tmpElmnt;
            spray_.deleteParticle(p1());
            p1 = tmpElmnt;
        }
        else 
        {
            ++p1;
        }
    } // outer loop

} // end


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
