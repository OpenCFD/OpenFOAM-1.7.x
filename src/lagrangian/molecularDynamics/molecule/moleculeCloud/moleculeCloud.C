/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

\*----------------------------------------------------------------------------*/

#include "moleculeCloud.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineParticleTypeNameAndDebug(molecule, 0);
    defineTemplateTypeNameAndDebug(Cloud<molecule>, 0);
};

Foam::scalar Foam::moleculeCloud::kb = 1.380650277e-23;

Foam::scalar Foam::moleculeCloud::elementaryCharge = 1.602176487e-19;

Foam::scalar Foam::moleculeCloud::vacuumPermittivity = 8.854187817e-12;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::moleculeCloud::buildConstProps()
{
    Info<< nl << "Reading moleculeProperties dictionary." << endl;

    const List<word>& idList(pot_.idList());

    constPropList_.setSize(idList.size());

    const List<word>& siteIdList(pot_.siteIdList());

    IOdictionary moleculePropertiesDict
    (
        IOobject
        (
            "moleculeProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    forAll(idList, i)
    {
        const word& id(idList[i]);

        const dictionary& molDict(moleculePropertiesDict.subDict(id));

        List<word> siteIdNames = molDict.lookup("siteIds");

        List<label> siteIds(siteIdNames.size());

        forAll(siteIdNames, sI)
        {
            const word& siteId = siteIdNames[sI];

            siteIds[sI] = findIndex(siteIdList, siteId);

            if (siteIds[sI] == -1)
            {
                FatalErrorIn("moleculeCloud.C") << nl
                    << siteId << " site not found."
                    << nl << abort(FatalError);
            }
        }

        molecule::constantProperties& constProp = constPropList_[i];

        constProp = molecule::constantProperties(molDict);

        constProp.siteIds() = siteIds;
    }
}


void Foam::moleculeCloud::setSiteSizesAndPositions()
{
    iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        const molecule::constantProperties& cP = constProps(mol().id());

        mol().setSiteSizes(cP.nSites());

        mol().setSitePositions(cP);
    }
}


void Foam::moleculeCloud::buildCellOccupancy()
{
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    iterator mol(this->begin());

    for
    (
        mol = this->begin();
        mol != this->end();
        ++mol
    )
    {
        cellOccupancy_[mol().cell()].append(&mol());
    }

    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].shrink();
    }

    il_.ril().referMolecules(cellOccupancy());
}


void Foam::moleculeCloud::calculatePairForce()
{
    iterator mol(this->begin());

    {
        // Real-Real interactions

        molecule* molI = &mol();

        molecule* molJ = &mol();

        const directInteractionList& dil(il_.dil());

        forAll(dil, d)
        {
            forAll(cellOccupancy_[d],cellIMols)
            {
                molI = cellOccupancy_[d][cellIMols];

                forAll(dil[d], interactingCells)
                {
                    List< molecule* > cellJ =
                    cellOccupancy_[dil[d][interactingCells]];

                    forAll(cellJ, cellJMols)
                    {
                        molJ = cellJ[cellJMols];

                        evaluatePair(molI, molJ);
                    }
                }

                forAll(cellOccupancy_[d],cellIOtherMols)
                {
                    molJ = cellOccupancy_[d][cellIOtherMols];

                    if (molJ > molI)
                    {
                        evaluatePair(molI, molJ);
                    }
                }
            }
        }
    }

    {
        // Real-Referred interactions

        molecule* molK = &mol();

        referredCellList& ril(il_.ril());

        forAll(ril, r)
        {
            const List<label>& realCells = ril[r].realCellsForInteraction();

            forAll(ril[r], refMols)
            {
                referredMolecule* molL(&(ril[r][refMols]));

                forAll(realCells, rC)
                {
                    List<molecule*> cellK = cellOccupancy_[realCells[rC]];

                    forAll(cellK, cellKMols)
                    {
                        molK = cellK[cellKMols];

                        evaluatePair(molK, molL);
                    }
                }
            }
        }
    }
}


void Foam::moleculeCloud::calculateTetherForce()
{
    const tetherPotentialList& tetherPot(pot_.tetherPotentials());

    iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        if (mol().tethered())
        {
            vector rIT = mol().position() - mol().specialPosition();

            label idI = mol().id();

            scalar massI = constProps(idI).mass();

            vector fIT = tetherPot.force(idI, rIT);

            mol().a() += fIT/massI;

            mol().potentialEnergy() += tetherPot.energy(idI, rIT);

            // What to do here?
            // mol().rf() += rIT*fIT;
        }
    }
}


void Foam::moleculeCloud::calculateExternalForce()
{
    iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        mol().a() += pot_.gravity();
    }
}


void Foam::moleculeCloud::removeHighEnergyOverlaps()
{
    Info<< nl << "Removing high energy overlaps, limit = "
        << pot_.potentialEnergyLimit()
        << nl << "Removal order:";

    forAll(pot_.removalOrder(), rO)
    {
        Info<< ' ' << pot_.idList()[pot_.removalOrder()[rO]];
    }

    Info<< nl ;

    label initialSize = this->size();

    buildCellOccupancy();

    // Real-Real interaction
    iterator mol(this->begin());

    {
        mol = this->begin();

        molecule* molI = &mol();

        molecule* molJ = &mol();

        DynamicList<molecule*> molsToDelete;

        const directInteractionList& dil(il_.dil());

        forAll(dil, d)
        {
            forAll(cellOccupancy_[d],cellIMols)
            {
                molI = cellOccupancy_[d][cellIMols];

                forAll(dil[d], interactingCells)
                {
                    List< molecule* > cellJ =
                    cellOccupancy_[dil[d][interactingCells]];

                    forAll(cellJ, cellJMols)
                    {
                        molJ = cellJ[cellJMols];

                        if (evaluatePotentialLimit(molI, molJ))
                        {
                            label idI = molI->id();

                            label idJ = molJ->id();

                            if
                            (
                                idI == idJ
                             || findIndex(pot_.removalOrder(), idJ)
                                < findIndex(pot_.removalOrder(), idI)
                            )
                            {
                                if (findIndex(molsToDelete, molJ) == -1)
                                {
                                    molsToDelete.append(molJ);
                                }
                            }
                            else if (findIndex(molsToDelete, molI) == -1)
                            {
                                molsToDelete.append(molI);
                            }
                        }
                    }
                }
            }

            forAll(cellOccupancy_[d],cellIOtherMols)
            {
                molJ = cellOccupancy_[d][cellIOtherMols];

                if (molJ > molI)
                {
                    if (evaluatePotentialLimit(molI, molJ))
                    {
                        label idI = molI->id();

                        label idJ = molJ->id();

                        if
                        (
                            idI == idJ
                         || findIndex(pot_.removalOrder(), idJ)
                            < findIndex(pot_.removalOrder(), idI)
                        )
                        {
                            if (findIndex(molsToDelete, molJ) == -1)
                            {
                                molsToDelete.append(molJ);
                            }
                        }
                        else if (findIndex(molsToDelete, molI) == -1)
                        {
                            molsToDelete.append(molI);
                        }
                    }
                }
            }
        }

        forAll (molsToDelete, mTD)
        {
            deleteParticle(*(molsToDelete[mTD]));
        }
    }

    buildCellOccupancy();

    // Real-Referred interaction

    {
        molecule* molK = &mol();

        DynamicList<molecule*> molsToDelete;

        referredCellList& ril(il_.ril());

        forAll(ril, r)
        {
            referredCell& refCell = ril[r];

            forAll(refCell, refMols)
            {
                referredMolecule* molL = &(refCell[refMols]);

                List <label> realCells = refCell.realCellsForInteraction();

                forAll(realCells, rC)
                {
                    label cellK = realCells[rC];

                    List<molecule*> cellKMols = cellOccupancy_[cellK];

                    forAll(cellKMols, cKM)
                    {
                        molK = cellKMols[cKM];

                        if (evaluatePotentialLimit(molK, molL))
                        {
                            label idK = molK->id();

                            label idL = molL->id();

                            if
                            (
                                findIndex(pot_.removalOrder(), idK)
                                    < findIndex(pot_.removalOrder(), idL)
                            )
                            {
                                if (findIndex(molsToDelete, molK) == -1)
                                {
                                    molsToDelete.append(molK);
                                }
                            }

                            else if
                            (
                                findIndex(pot_.removalOrder(), idK)
                                    == findIndex(pot_.removalOrder(), idL)
                            )
                            {
                                if
                                (
                                    Pstream::myProcNo() == refCell.sourceProc()
                                 && cellK <= refCell.sourceCell()
                                )
                                {
                                    if (findIndex(molsToDelete, molK) == -1)
                                    {
                                        molsToDelete.append(molK);
                                    }
                                }

                                else if
                                (
                                    Pstream::myProcNo() < refCell.sourceProc()
                                )
                                {
                                    if (findIndex(molsToDelete, molK) == -1)
                                    {
                                        molsToDelete.append(molK);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        forAll (molsToDelete, mTD)
        {
            deleteParticle(*(molsToDelete[mTD]));
        }
    }

    buildCellOccupancy();

    label molsRemoved = initialSize - this->size();

    if (Pstream::parRun())
    {
        reduce(molsRemoved, sumOp<label>());
    }

    Info<< tab << molsRemoved << " molecules removed" << endl;
}


void Foam::moleculeCloud::initialiseMolecules
(
    const IOdictionary& mdInitialiseDict
)
{
    Info<< nl
        << "Initialising molecules in each zone specified in mdInitialiseDict."
        << endl;

    const cellZoneMesh& cellZones = mesh_.cellZones();

    if (!cellZones.size())
    {
        FatalErrorIn("void Foam::moleculeCloud::initialiseMolecules")
            << "No cellZones found in the mesh."
            << abort(FatalError);
    }

    forAll (cellZones, z)
    {
        const cellZone& zone(cellZones[z]);

        if (zone.size())
        {
            if (!mdInitialiseDict.found(zone.name()))
            {
                Info<< "No specification subDictionary for zone "
                    << zone.name() << " found, skipping." << endl;
            }
            else
            {
                const dictionary& zoneDict =
                    mdInitialiseDict.subDict(zone.name());

                const scalar temperature
                (
                    readScalar(zoneDict.lookup("temperature"))
                );

                const vector bulkVelocity(zoneDict.lookup("bulkVelocity"));

                List<word> latticeIds
                (
                    zoneDict.lookup("latticeIds")
                );

                List<vector> latticePositions
                (
                    zoneDict.lookup("latticePositions")
                );

                if(latticeIds.size() != latticePositions.size())
                {
                    FatalErrorIn("Foam::moleculeCloud::initialiseMolecules")
                        << "latticeIds and latticePositions must be the same "
                        << " size." << nl
                        << abort(FatalError);
                }

                diagTensor latticeCellShape
                (
                    zoneDict.lookup("latticeCellShape")
                );

                scalar latticeCellScale = 0.0;

                if (zoneDict.found("numberDensity"))
                {
                    scalar numberDensity = readScalar
                    (
                        zoneDict.lookup("numberDensity")
                    );

                    if (numberDensity < VSMALL)
                    {
                        WarningIn("moleculeCloud::initialiseMolecules")
                            << "numberDensity too small, not filling zone "
                            << zone.name() << endl;

                        continue;
                    }

                    latticeCellScale = pow
                    (
                        latticeIds.size()/(det(latticeCellShape)*numberDensity),
                        (1.0/3.0)
                    );
                }
                else if (zoneDict.found("massDensity"))
                {
                    scalar unitCellMass = 0.0;

                    forAll(latticeIds, i)
                    {
                        label id = findIndex(pot_.idList(), latticeIds[i]);

                        const molecule::constantProperties& cP(constProps(id));

                        unitCellMass += cP.mass();
                    }

                    scalar massDensity = readScalar
                    (
                        zoneDict.lookup("massDensity")
                    );

                    if (massDensity < VSMALL)
                    {
                        WarningIn("moleculeCloud::initialiseMolecules")
                            << "massDensity too small, not filling zone "
                            << zone.name() << endl;

                        continue;
                    }


                    latticeCellScale = pow
                    (
                        unitCellMass/(det(latticeCellShape)*massDensity),
                        (1.0/3.0)
                    );
                }
                else
                {
                    FatalErrorIn("Foam::moleculeCloud::initialiseMolecules")
                        << "massDensity or numberDensity not specified " << nl
                        << abort(FatalError);
                }

                latticeCellShape *= latticeCellScale;

                vector anchor(zoneDict.lookup("anchor"));

                bool tethered = false;

                if (zoneDict.found("tetherSiteIds"))
                {
                    tethered = bool
                    (
                        List<word>(zoneDict.lookup("tetherSiteIds")).size()
                    );
                }

                const vector orientationAngles
                (
                    zoneDict.lookup("orientationAngles")
                );

                scalar phi
                (
                    orientationAngles.x()*mathematicalConstant::pi/180.0
                );

                scalar theta
                (
                    orientationAngles.y()*mathematicalConstant::pi/180.0
                );

                scalar psi
                (
                    orientationAngles.z()*mathematicalConstant::pi/180.0
                );

                const tensor R
                (
                    cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
                    cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
                    sin(psi)*sin(theta),
                    - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                    - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                    cos(psi)*sin(theta),
                    sin(theta)*sin(phi),
                    - sin(theta)*cos(phi),
                    cos(theta)
                );

                // Find the optimal anchor position.  Finding the approximate
                // mid-point of the zone of cells and snapping to the nearest
                // lattice location.

                vector zoneMin = VGREAT*vector::one;

                vector zoneMax = -VGREAT*vector::one;

                forAll (zone, cell)
                {
                    const point cellCentre = mesh_.cellCentres()[zone[cell]];

                    if (cellCentre.x() > zoneMax.x())
                    {
                        zoneMax.x() = cellCentre.x();
                    }
                    if (cellCentre.x() < zoneMin.x())
                    {
                        zoneMin.x() = cellCentre.x();
                    }
                    if (cellCentre.y() > zoneMax.y())
                    {
                        zoneMax.y() = cellCentre.y();
                    }
                    if (cellCentre.y() < zoneMin.y())
                    {
                        zoneMin.y() = cellCentre.y();
                    }
                    if (cellCentre.z() > zoneMax.z())
                    {
                        zoneMax.z() = cellCentre.z();
                    }
                    if (cellCentre.z() < zoneMin.z())
                    {
                        zoneMin.z() = cellCentre.z();
                    }
                }

                point zoneMid = 0.5*(zoneMin + zoneMax);

                point latticeMid = inv(latticeCellShape) & (R.T()
                & (zoneMid - anchor));

                point latticeAnchor
                (
                    label(latticeMid.x() + 0.5*sign(latticeMid.x())),
                    label(latticeMid.y() + 0.5*sign(latticeMid.y())),
                    label(latticeMid.z() + 0.5*sign(latticeMid.z()))
                );

                anchor += (R & (latticeCellShape & latticeAnchor));

                // Continue trying to place molecules as long as at
                // least one molecule is placed in each iteration.
                // The "|| totalZoneMols == 0" condition means that the
                // algorithm will continue if the origin is outside the
                // zone.

                label n = 0;

                label totalZoneMols = 0;

                label molsPlacedThisIteration = 0;

                while
                (
                    molsPlacedThisIteration != 0
                 || totalZoneMols == 0
                )
                {
                    label sizeBeforeIteration = this->size();

                    bool partOfLayerInBounds = false;

                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    // start of placement of molecules
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    if (n == 0)
                    {
                        // Special treatment is required for the first position,
                        // i.e. iteration zero.

                        labelVector unitCellLatticePosition(0,0,0);

                        forAll(latticePositions, p)
                        {
                            label id = findIndex(pot_.idList(), latticeIds[p]);

                            const vector& latticePosition =
                            vector
                            (
                                unitCellLatticePosition.x(),
                                unitCellLatticePosition.y(),
                                unitCellLatticePosition.z()
                            )
                            + latticePositions[p];

                            point globalPosition =
                            anchor + (R & (latticeCellShape & latticePosition));

                            partOfLayerInBounds =
                            mesh_.bounds().contains(globalPosition);

                            label cell = mesh_.findCell(globalPosition);

                            if (findIndex(zone, cell) != -1)
                            {
                                createMolecule
                                (
                                    globalPosition,
                                    cell,
                                    id,
                                    tethered,
                                    temperature,
                                    bulkVelocity
                                );
                            }
                        }
                    }
                    else
                    {
                        // Place top and bottom caps.

                        labelVector unitCellLatticePosition(0,0,0);

                        for
                        (
                            unitCellLatticePosition.z() = -n;
                            unitCellLatticePosition.z() <= n;
                            unitCellLatticePosition.z() += 2*n
                        )
                        {
                            for
                            (
                                unitCellLatticePosition.y() = -n;
                                unitCellLatticePosition.y() <= n;
                                unitCellLatticePosition.y()++
                            )
                            {
                                for
                                (
                                    unitCellLatticePosition.x() = -n;
                                    unitCellLatticePosition.x() <= n;
                                    unitCellLatticePosition.x()++
                                )
                                {
                                    forAll(latticePositions, p)
                                    {
                                        label id = findIndex
                                        (
                                            pot_.idList(),
                                            latticeIds[p]
                                        );

                                        const vector& latticePosition =
                                        vector
                                        (
                                            unitCellLatticePosition.x(),
                                            unitCellLatticePosition.y(),
                                            unitCellLatticePosition.z()
                                        )
                                        + latticePositions[p];

                                        point globalPosition = anchor
                                        + (R
                                        &(latticeCellShape & latticePosition));

                                        partOfLayerInBounds =
                                        mesh_.bounds().contains(globalPosition);

                                        label cell = mesh_.findCell
                                        (
                                            globalPosition
                                        );

                                        if (findIndex(zone, cell) != -1)
                                        {
                                            createMolecule
                                            (
                                                globalPosition,
                                                cell,
                                                id,
                                                tethered,
                                                temperature,
                                                bulkVelocity
                                            );
                                        }
                                    }
                                }
                            }
                        }

                        for
                        (
                            unitCellLatticePosition.z() = -(n-1);
                            unitCellLatticePosition.z() <= (n-1);
                            unitCellLatticePosition.z()++
                        )
                        {
                            for (label iR = 0; iR <= 2*n -1; iR++)
                            {
                                unitCellLatticePosition.x() = n;

                                unitCellLatticePosition.y() = -n + (iR + 1);

                                for (label iK = 0; iK < 4; iK++)
                                {
                                    forAll(latticePositions, p)
                                    {
                                        label id = findIndex
                                        (
                                            pot_.idList(),
                                            latticeIds[p]
                                        );

                                        const vector& latticePosition =
                                        vector
                                        (
                                            unitCellLatticePosition.x(),
                                            unitCellLatticePosition.y(),
                                            unitCellLatticePosition.z()
                                        )
                                        + latticePositions[p];

                                        point globalPosition = anchor
                                        + (R
                                        &(latticeCellShape & latticePosition));

                                        partOfLayerInBounds =
                                        mesh_.bounds().contains(globalPosition);

                                        label cell = mesh_.findCell
                                        (
                                            globalPosition
                                        );

                                        if (findIndex(zone, cell) != -1)
                                        {
                                            createMolecule
                                            (
                                                globalPosition,
                                                cell,
                                                id,
                                                tethered,
                                                temperature,
                                                bulkVelocity
                                            );
                                        }
                                    }

                                    unitCellLatticePosition =
                                    labelVector
                                    (
                                        - unitCellLatticePosition.y(),
                                        unitCellLatticePosition.x(),
                                        unitCellLatticePosition.z()
                                    );
                                }
                            }
                        }
                    }

                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    // end of placement of molecules
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    if
                    (
                        totalZoneMols == 0
                     && !partOfLayerInBounds
                    )
                    {
                        WarningIn("Foam::moleculeCloud::initialiseMolecules()")
                            << "A whole layer of unit cells was placed "
                            << "outside the bounds of the mesh, but no "
                            << "molecules have been placed in zone '"
                            << zone.name()
                            << "'.  This is likely to be because the zone "
                            << "has few cells ("
                            << zone.size()
                            << " in this case) and no lattice position "
                            << "fell inside them.  "
                            << "Aborting filling this zone."
                            << endl;

                        break;
                    }

                    molsPlacedThisIteration =
                    this->size() - sizeBeforeIteration;

                    totalZoneMols += molsPlacedThisIteration;

                    n++;
                }
            }
        }
    }
}


void Foam::moleculeCloud::createMolecule
(
    const point& position,
    label cell,
    label id,
    bool tethered,
    scalar temperature,
    const vector& bulkVelocity
)
{
    const Cloud<molecule>& cloud = *this;

    if (cell == -1)
    {
        cell = mesh_.findCell(position);
    }

    if (cell == -1)
    {
        FatalErrorIn("Foam::moleculeCloud::createMolecule")
            << "Position specified does not correspond to a mesh cell." << nl
            << abort(FatalError);
    }

    point specialPosition(vector::zero);

    label special = 0;

    if (tethered)
    {
        specialPosition = position;

        special = molecule::SPECIAL_TETHERED;
    }

    const molecule::constantProperties& cP(constProps(id));

    vector v = equipartitionLinearVelocity(temperature, cP.mass());

    v += bulkVelocity;

    vector pi = vector::zero;

    tensor Q = I;

    if (!cP.pointMolecule())
    {
        pi = equipartitionAngularMomentum(temperature, cP);

        scalar phi(rndGen_.scalar01()*mathematicalConstant::twoPi);

        scalar theta(rndGen_.scalar01()*mathematicalConstant::twoPi);

        scalar psi(rndGen_.scalar01()*mathematicalConstant::twoPi);

        Q = tensor
        (
            cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
            cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
            sin(psi)*sin(theta),
            - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
            - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
            cos(psi)*sin(theta),
            sin(theta)*sin(phi),
            - sin(theta)*cos(phi),
            cos(theta)
        );
    }

    addParticle
    (
        new molecule
        (
            cloud,
            position,
            cell,
            Q,
            v,
            vector::zero,
            pi,
            vector::zero,
            specialPosition,
            constProps(id),
            special,
            id
        )
    );
}


Foam::label Foam::moleculeCloud::nSites() const
{
    label n = 0;

    const_iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        n += constProps(mol().id()).nSites();
    }

    return n;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::moleculeCloud::moleculeCloud
(
    const polyMesh& mesh,
    const potential& pot,
    bool readFields
)
:
    Cloud<molecule>(mesh, "moleculeCloud", false),
    mesh_(mesh),
    pot_(pot),
    cellOccupancy_(mesh_.nCells()),
    il_(mesh_, pot_.pairPotentials().rCutMaxSqr(), false),
    constPropList_(),
    rndGen_(clock::getTime())
{
    if (readFields)
    {
        molecule::readFields(*this);
    }

    buildConstProps();

    setSiteSizesAndPositions();

    removeHighEnergyOverlaps();

    calculateForce();
}


Foam::moleculeCloud::moleculeCloud
(
    const polyMesh& mesh,
    const potential& pot,
    const IOdictionary& mdInitialiseDict,
    bool readFields
)
    :
    Cloud<molecule>(mesh, "moleculeCloud", false),
    mesh_(mesh),
    pot_(pot),
    il_(mesh_),
    constPropList_(),
    rndGen_(clock::getTime())
{
    if (readFields)
    {
        molecule::readFields(*this);
    }

    clear();

    buildConstProps();

    initialiseMolecules(mdInitialiseDict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::moleculeCloud::evolve()
{
    molecule::trackData td0(*this, 0);
    Cloud<molecule>::move(td0);

    molecule::trackData td1(*this, 1);
    Cloud<molecule>::move(td1);

    molecule::trackData td2(*this, 2);
    Cloud<molecule>::move(td2);

    calculateForce();

    molecule::trackData td3(*this, 3);
    Cloud<molecule>::move(td3);
}


void Foam::moleculeCloud::calculateForce()
{
    buildCellOccupancy();

    iterator mol(this->begin());

    // Set accumulated quantities to zero
    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        mol().siteForces() = vector::zero;

        mol().potentialEnergy() = 0.0;

        mol().rf() = tensor::zero;
    }

    calculatePairForce();

    calculateTetherForce();

    calculateExternalForce();
}


void Foam::moleculeCloud::applyConstraintsAndThermostats
(
    const scalar targetTemperature,
    const scalar measuredTemperature
)
{
    scalar temperatureCorrectionFactor =
        sqrt(targetTemperature/measuredTemperature);

    Info<< "----------------------------------------" << nl
        << "Temperature equilibration" << nl
        << "Target temperature = "
        << targetTemperature << nl
        << "Measured temperature = "
        << measuredTemperature << nl
        << "Temperature correction factor = "
        << temperatureCorrectionFactor << nl
        << "----------------------------------------"
        << endl;

    iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        mol().v() *= temperatureCorrectionFactor;

        mol().pi() *= temperatureCorrectionFactor;
    }
}


void Foam::moleculeCloud::writeXYZ(const fileName& fName) const
{
    OFstream str(fName);

    str << nSites() << nl << "moleculeCloud site positions in angstroms" << nl;

    const_iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        const molecule::constantProperties& cP = constProps(mol().id());

        forAll(mol().sitePositions(), i)
        {
            const point& sP = mol().sitePositions()[i];

            str << pot_.siteIdList()[cP.siteIds()[i]]
                << ' ' << sP.x()*1e10
                << ' ' << sP.y()*1e10
                << ' ' << sP.z()*1e10
                << nl;
        }
    }
}


// ************************************************************************* //
