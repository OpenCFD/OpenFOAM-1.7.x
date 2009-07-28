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

#include "processorFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(processorFvPatch, 0);
addToRunTimeSelectionTable(fvPatch, processorFvPatch, polyPatch);


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void processorFvPatch::makeWeights(scalarField& w) const
{
    if (Pstream::parRun())
    {
        // The face normals point in the opposite direction on the other side
        scalarField neighbFaceCentresCn
        (
            (
                procPolyPatch_.neighbFaceAreas()
               /(mag(procPolyPatch_.neighbFaceAreas()) + VSMALL)
            )
          & (
              procPolyPatch_.neighbFaceCentres()
            - procPolyPatch_.neighbFaceCellCentres())
        );

        w = neighbFaceCentresCn/((nf()&fvPatch::delta()) + neighbFaceCentresCn);
    }
    else
    {
        w = 1.0;
    }
}


void processorFvPatch::makeDeltaCoeffs(scalarField& dc) const
{
    if (Pstream::parRun())
    {
        dc = (1.0 - weights())/(nf() & fvPatch::delta());
    }
    else
    {
        dc = 1.0/(nf() & fvPatch::delta());
    }
}


tmp<vectorField> processorFvPatch::delta() const
{
    if (Pstream::parRun())
    {
        // To the transformation if necessary
        if (parallel())
        {
            return
                fvPatch::delta()
              - (
                    procPolyPatch_.neighbFaceCentres()
                  - procPolyPatch_.neighbFaceCellCentres()
                );
        }
        else
        {
            return
                fvPatch::delta()
              - transform
                (
                    forwardT(),
                    (
                        procPolyPatch_.neighbFaceCentres()
                      - procPolyPatch_.neighbFaceCellCentres()
                    )
                );
        }
    }
    else
    {
        return fvPatch::delta();
    }
}


tmp<labelField> processorFvPatch::interfaceInternalField
(
    const unallocLabelList& internalData
) const
{
    return patchInternalField(internalData);
}


void processorFvPatch::initTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& interfaceData
) const
{
    send(commsType, interfaceData);
}


tmp<labelField> processorFvPatch::transfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList&
) const
{
    return receive<label>(commsType, this->size());
}


void processorFvPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList& iF
) const
{
    send(commsType, patchInternalField(iF)());
}


tmp<labelField> processorFvPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const unallocLabelList&
) const
{
    return receive<label>(commsType, this->size());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
