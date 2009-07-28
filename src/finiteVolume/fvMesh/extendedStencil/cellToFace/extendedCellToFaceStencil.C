/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "extendedCellToFaceStencil.H"
#include "globalIndex.H"
#include "syncTools.H"
#include "SortableList.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

defineTypeNameAndDebug(Foam::extendedCellToFaceStencil, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::extendedCellToFaceStencil::writeStencilStats
(
    Ostream& os,
    const labelListList& stencil,
    const mapDistribute& map
)
{
    label sumSize = 0;
    label nSum = 0;
    label minSize = labelMax;
    label maxSize = labelMin;

    forAll(stencil, i)
    {
        const labelList& sCells = stencil[i];

        if (sCells.size() > 0)
        {
            sumSize += sCells.size();
            nSum++;
            minSize = min(minSize, sCells.size());
            maxSize = max(maxSize, sCells.size());
        }
    }
    reduce(sumSize, sumOp<label>());
    reduce(nSum, sumOp<label>());

    reduce(minSize, minOp<label>());
    reduce(maxSize, maxOp<label>());

    os  << "Stencil size :" << nl
        << "    average : " << scalar(sumSize)/nSum << nl
        << "    min     : " << minSize << nl
        << "    max     : " << maxSize << nl
        << endl;

    // Sum all sent data
    label nSent = 0;
    label nLocal = 0;
    forAll(map.subMap(), procI)
    {
        if (procI != Pstream::myProcNo())
        {
            nSent += map.subMap()[procI].size();
        }
        else
        {
            nLocal += map.subMap()[procI].size();
        }
    }

    os  << "Local data size : " << returnReduce(nLocal, sumOp<label>()) << nl
        << "Sent data size  : " << returnReduce(nSent, sumOp<label>()) << nl
        << endl;
}


Foam::autoPtr<Foam::mapDistribute>
Foam::extendedCellToFaceStencil::calcDistributeMap
(
    const polyMesh& mesh,
    const globalIndex& globalNumbering,
    labelListList& faceStencil
)
{
    // Convert stencil to schedule
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // We now know what information we need from other processors. This needs
    // to be converted into what information I need to send as well
    // (mapDistribute)


    // 1. Construct per processor compact addressing of the global cells
    //    needed. The ones from the local processor are not included since
    //    these are always all needed.
    List<Map<label> > globalToProc(Pstream::nProcs());
    {
        const labelList& procPatchMap = mesh.globalData().procPatchMap();
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // Presize with (as estimate) size of patch to neighbour.
        forAll(procPatchMap, procI)
        {
            if (procPatchMap[procI] != -1)
            {
                globalToProc[procI].resize
                (
                    patches[procPatchMap[procI]].size()
                );
            }
        }

        // Collect all (non-local) globalcells/faces needed.
        forAll(faceStencil, faceI)
        {
            const labelList& stencilCells = faceStencil[faceI];

            forAll(stencilCells, i)
            {
                label globalCellI = stencilCells[i];
                label procI = globalNumbering.whichProcID(stencilCells[i]);

                if (procI != Pstream::myProcNo())
                {
                    label nCompact = globalToProc[procI].size();
                    globalToProc[procI].insert(globalCellI, nCompact);
                }
            }
        }
        // Sort global cells needed (not really necessary)
        forAll(globalToProc, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                Map<label>& globalMap = globalToProc[procI];

                SortableList<label> sorted(globalMap.toc().xfer());

                forAll(sorted, i)
                {
                    Map<label>::iterator iter = globalMap.find(sorted[i]);
                    iter() = i;
                }
            }
        }


        //forAll(globalToProc, procI)
        //{
        //    Pout<< "From processor:" << procI << " want cells/faces:" << endl;
        //    forAllConstIter(Map<label>, globalToProc[procI], iter)
        //    {
        //        Pout<< "    global:" << iter.key()
        //            << " local:" << globalNumbering.toLocal(procI, iter.key())
        //            << endl;
        //    }
        //    Pout<< endl;
        //}
    }


    // 2. The overall compact addressing is
    // - myProcNo data first (uncompacted)
    // - all other processors consecutively

    labelList compactStart(Pstream::nProcs());
    compactStart[Pstream::myProcNo()] = 0;
    label nCompact = globalNumbering.localSize();
    forAll(compactStart, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            compactStart[procI] = nCompact;
            nCompact += globalToProc[procI].size();
        }
    }


    // 3. Find out what to receive/send in compact addressing.
    labelListList recvCompact(Pstream::nProcs());
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            labelList wantedGlobals(globalToProc[procI].size());
            recvCompact[procI].setSize(globalToProc[procI].size());

            label i = 0;
            forAllConstIter(Map<label>, globalToProc[procI], iter)
            {
                wantedGlobals[i] = iter.key();
                recvCompact[procI][i] = compactStart[procI]+iter();
                i++;
            }

            // Send the global cell numbers I need from procI
            OPstream str(Pstream::blocking, procI);
            str << wantedGlobals;
        }
        else
        {
            recvCompact[procI] =
                compactStart[procI]
              + identity(globalNumbering.localSize());
        }
    }
    labelListList sendCompact(Pstream::nProcs());
    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        if (procI != Pstream::myProcNo())
        {
            // See what neighbour wants to receive (= what I need to send)

            IPstream str(Pstream::blocking, procI);
            labelList globalCells(str);

            labelList& procCompact = sendCompact[procI];
            procCompact.setSize(globalCells.size());

            // Convert from globalCells (all on my processor!) into compact
            // addressing
            forAll(globalCells, i)
            {
                label cellI = globalNumbering.toLocal(globalCells[i]);
                procCompact[i] = compactStart[Pstream::myProcNo()]+cellI;
            }
        }
        else
        {
            sendCompact[procI] = recvCompact[procI];
        }
    }

    // Convert stencil to compact numbering
    forAll(faceStencil, faceI)
    {
        labelList& stencilCells = faceStencil[faceI];

        forAll(stencilCells, i)
        {
            label globalCellI = stencilCells[i];
            label procI = globalNumbering.whichProcID(globalCellI);
            if (procI != Pstream::myProcNo())
            {
                label localCompact = globalToProc[procI][globalCellI];
                stencilCells[i] = compactStart[procI]+localCompact;
            }
            else
            {
                label localCompact = globalNumbering.toLocal(globalCellI);
                stencilCells[i] = compactStart[procI]+localCompact;
            }

        }
    }


    // Constuct map for distribution of compact data.
    autoPtr<mapDistribute> mapPtr
    (
        new mapDistribute
        (
            nCompact,
            sendCompact,
            recvCompact,
            true            // reuse send/recv maps.
        )
    );

    return mapPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::extendedCellToFaceStencil::extendedCellToFaceStencil(const polyMesh& mesh)
:
    mesh_(mesh)
{
    // Check for transformation - not supported.
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        if (isA<coupledPolyPatch>(patches[patchI]))
        {
            const coupledPolyPatch& cpp =
                refCast<const coupledPolyPatch>(patches[patchI]);

            if (!cpp.parallel() || cpp.separated())
            {
                FatalErrorIn
                (
                    "extendedCellToFaceStencil::extendedCellToFaceStencil"
                    "(const polyMesh&)"
                )   << "Coupled patches with transformations not supported."
                    << endl
                    << "Problematic patch " << cpp.name() << exit(FatalError);
            }
        }
    }
}


// ************************************************************************* //
