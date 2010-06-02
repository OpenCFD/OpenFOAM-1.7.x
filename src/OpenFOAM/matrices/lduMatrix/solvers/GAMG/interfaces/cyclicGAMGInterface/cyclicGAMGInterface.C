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

#include "cyclicGAMGInterface.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicGAMGInterface, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterface,
        cyclicGAMGInterface,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicGAMGInterface::cyclicGAMGInterface
(
    const lduInterface& fineInterface,
    const labelField& localRestrictAddressing,
    const labelField& neighbourRestrictAddressing
)
:
    GAMGInterface
    (
        fineInterface,
        localRestrictAddressing,
        neighbourRestrictAddressing
    ),
    fineCyclicInterface_(refCast<const cyclicLduInterface>(fineInterface))
{
    // Make a lookup table of entries for owner/neighbour
    HashTable<SLList<label>, label, Hash<label> > neighboursTable
    (
        localRestrictAddressing.size()
    );

    // Table of face-sets to be agglomerated
    HashTable<SLList<SLList<label> >, label, Hash<label> > faceFaceTable
    (
        localRestrictAddressing.size()
    );

    label nCoarseFaces = 0;

    label sizeBy2 = localRestrictAddressing.size()/2;

    for (label ffi=0; ffi<sizeBy2; ffi++)
    {
        label curMaster = localRestrictAddressing[ffi];
        label curSlave = localRestrictAddressing[ffi + sizeBy2];

        // Look for the master cell.  If it has already got a face,
        // add the coefficient to the face.  If not, create a new
        // face.
        if (neighboursTable.found(curMaster))
        {
            // Check all current neighbours to see if the current
            // slave already exists.  If so, add the coefficient.

            SLList<label>& curNbrs = neighboursTable.find(curMaster)();

            SLList<SLList<label> >& curFaceFaces =
                faceFaceTable.find(curMaster)();

            bool nbrFound = false;

            SLList<label>::iterator nbrsIter = curNbrs.begin();

            SLList<SLList<label> >::iterator faceFacesIter =
                curFaceFaces.begin();

            for
            (
                ;
                nbrsIter != curNbrs.end(), faceFacesIter != curFaceFaces.end();
                ++nbrsIter, ++faceFacesIter
            )
            {
                if (nbrsIter() == curSlave)
                {
                    nbrFound = true;
                    faceFacesIter().append(ffi);
                    break;
                }
            }

            if (!nbrFound)
            {
                curNbrs.append(curSlave);
                curFaceFaces.append(ffi);

                // New coarse face created
                nCoarseFaces++;
            }
        }
        else
        {
            // This master has got no neighbours yet.  Add a neighbour
            // and a coefficient, thus creating a new face
            neighboursTable.insert(curMaster, SLList<label>(curSlave));
            faceFaceTable.insert(curMaster, SLList<SLList<label> >(ffi));

            // New coarse face created
            nCoarseFaces++;
        }
    } // end for all fine faces


    faceCells_.setSize(2*nCoarseFaces, -1);
    faceRestrictAddressing_.setSize(localRestrictAddressing.size(), -1);

    labelList contents = neighboursTable.toc();

    // Reset face counter for re-use
    nCoarseFaces = 0;

    // On master side, the owner addressing is stored in table of contents
    forAll (contents, masterI)
    {
        SLList<label>& curNbrs = neighboursTable.find(contents[masterI])();

        SLList<SLList<label> >& curFaceFaces =
            faceFaceTable.find(contents[masterI])();

        SLList<label>::iterator nbrsIter = curNbrs.begin();
        SLList<SLList<label> >::iterator faceFacesIter = curFaceFaces.begin();

        for
        (
            ;
            nbrsIter != curNbrs.end(), faceFacesIter != curFaceFaces.end();
            ++nbrsIter, ++faceFacesIter
        )
        {
            faceCells_[nCoarseFaces] = contents[masterI];

            for
            (
                SLList<label>::iterator facesIter = faceFacesIter().begin();
                facesIter != faceFacesIter().end();
                ++facesIter
            )
            {
                faceRestrictAddressing_[facesIter()] = nCoarseFaces;
            }

            nCoarseFaces++;
        }
    }

    // On slave side, the owner addressing is stored in linked lists
    forAll (contents, masterI)
    {
        SLList<label>& curNbrs = neighboursTable.find(contents[masterI])();

        SLList<SLList<label> >& curFaceFaces =
            faceFaceTable.find(contents[masterI])();

        SLList<label>::iterator nbrsIter = curNbrs.begin();
        SLList<SLList<label> >::iterator faceFacesIter = curFaceFaces.begin();

        for
        (
            ;
            nbrsIter != curNbrs.end(), faceFacesIter != curFaceFaces.end();
            ++nbrsIter, ++faceFacesIter
        )
        {
            faceCells_[nCoarseFaces] = nbrsIter();

            for
            (
                SLList<label>::iterator facesIter = faceFacesIter().begin();
                facesIter != faceFacesIter().end();
                ++facesIter
            )
            {
                faceRestrictAddressing_[facesIter() + sizeBy2] = nCoarseFaces;
            }

            nCoarseFaces++;
        }
    }
}


// * * * * * * * * * * * * * * * * Desstructor * * * * * * * * * * * * * * * //

Foam::cyclicGAMGInterface::~cyclicGAMGInterface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::labelField> Foam::cyclicGAMGInterface::transfer
(
    const Pstream::commsTypes,
    const unallocLabelList& interfaceData
) const
{
    tmp<labelField> tpnf(new labelField(size()));
    labelField& pnf = tpnf();

    label sizeby2 = size()/2;

    for (label facei=0; facei<sizeby2; facei++)
    {
        pnf[facei] = interfaceData[facei + sizeby2];
        pnf[facei + sizeby2] = interfaceData[facei];
    }

    return tpnf;
}


Foam::tmp<Foam::labelField> Foam::cyclicGAMGInterface::internalFieldTransfer
(
    const Pstream::commsTypes,
    const unallocLabelList& iF
) const
{
    tmp<labelField> tpnf(new labelField(size()));
    labelField& pnf = tpnf();

    label sizeby2 = size()/2;

    for (label facei=0; facei<sizeby2; facei++)
    {
        pnf[facei] = iF[faceCells_[facei + sizeby2]];
        pnf[facei + sizeby2] = iF[faceCells_[facei]];
    }

    return tpnf;
}


// ************************************************************************* //
