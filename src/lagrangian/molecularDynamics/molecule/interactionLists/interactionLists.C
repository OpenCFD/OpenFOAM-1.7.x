/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "interactionLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::interactionLists::transTol = 1e-12;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interactionLists::buildCellReferralLists()
{
    Info<< nl << "Determining molecule referring schedule" << endl;

    const referredCellList& refIntL(ril());

    DynamicList<label> referralProcs;

    // Run through all referredCells to build list of interacting processors

    forAll(refIntL, rIL)
    {
        const referredCell& rC(refIntL[rIL]);

        if (findIndex(referralProcs, rC.sourceProc()) == -1)
        {
            referralProcs.append(rC.sourceProc());
        }
    }

    referralProcs.shrink();

//     Pout << "referralProcs: " << nl << referralProcs << endl;

    List<DynamicList<label> > cellSendingReferralLists(referralProcs.size());

    List<DynamicList<DynamicList<label> > >
        cellReceivingReferralLists(referralProcs.size());

    // Run through all referredCells again building up send and receive info

    forAll(refIntL, rIL)
    {
        const referredCell& rC(refIntL[rIL]);

        label rPI = findIndex(referralProcs, rC.sourceProc());

        DynamicList<DynamicList<label> >& rRL(cellReceivingReferralLists[rPI]);

        DynamicList<label>& sRL(cellSendingReferralLists[rPI]);

        label existingSource = findIndex(sRL, rC.sourceCell());

        // Check to see if this source cell has already been allocated to
        // come to this processor.  If not, add the source cell to the sending
        // list and add the current referred cell to the receiving list.

        // It shouldn't be possible for the sending and receiving lists to be
        // different lengths, because their append operations happen at the
        // same time.

        if (existingSource == -1)
        {
            sRL.append(rC.sourceCell());

            rRL.append
            (
                DynamicList<label> (labelList(1,rIL))
            );
        }
        else
        {
            rRL[existingSource].append(rIL);

            rRL[existingSource].shrink();
        }
    }

    forAll(referralProcs, rPI)
    {
        DynamicList<DynamicList<label> >& rRL(cellReceivingReferralLists[rPI]);

        DynamicList<label>& sRL(cellSendingReferralLists[rPI]);

        sRL.shrink();

        rRL.shrink();
    }

    // It is assumed that cell exchange is reciprocal, if proc A has cells to
    // send to proc B, then proc B must have some to send to proc A.

    cellReceivingReferralLists_.setSize(referralProcs.size());

    cellSendingReferralLists_.setSize(referralProcs.size());

    forAll(referralProcs, rPI)
    {
        DynamicList<DynamicList<label> >& rRL(cellReceivingReferralLists[rPI]);

        labelListList translLL(rRL.size());

        forAll(rRL, rRLI)
        {
            translLL[rRLI] = rRL[rRLI];
        }

        cellReceivingReferralLists_[rPI] = receivingReferralList
        (
            referralProcs[rPI],
            translLL
        );
    }

    // Send sendingReferralLists to each interacting processor.

    forAll(referralProcs, rPI)
    {

        DynamicList<label>& sRL(cellSendingReferralLists[rPI]);

        if (referralProcs[rPI] != Pstream::myProcNo())
        {
            if (Pstream::parRun())
            {
                OPstream toInteractingProc
                (
                    Pstream::blocking,
                    referralProcs[rPI]
                );

                toInteractingProc << sendingReferralList
                (
                    Pstream::myProcNo(),
                    sRL
                );
            }
        }
    }

    // Receive sendingReferralLists from each interacting processor.

    forAll(referralProcs, rPI)
    {
        if (referralProcs[rPI] != Pstream::myProcNo())
        {
            if (Pstream::parRun())
            {
                IPstream fromInteractingProc
                (
                    Pstream::blocking,
                    referralProcs[rPI]
                );

                fromInteractingProc >> cellSendingReferralLists_[rPI];
            }
        }
        else
        {
            cellSendingReferralLists_[rPI] = sendingReferralList
            (
                Pstream::myProcNo(),
                cellSendingReferralLists[rPI]
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::interactionLists::interactionLists
(
    const polyMesh& mesh,
    scalar rCutMaxSqr,
    bool pointPointListBuild
)
:
    mesh_(mesh),
    rCutMaxSqr_(rCutMaxSqr),
    dil_(*this, pointPointListBuild),
    ril_(*this, pointPointListBuild),
    cellSendingReferralLists_(),
    cellReceivingReferralLists_()
{
    buildCellReferralLists();
}


Foam::interactionLists::interactionLists(const polyMesh& mesh)
:
    mesh_(mesh),
    dil_(*this),
    ril_(*this)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interactionLists::~interactionLists()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::interactionLists::testPointPointDistance
(
    const label ptI,
    const label ptJ
) const
{
    return (magSqr(mesh_.points()[ptI] - mesh_.points()[ptJ]) <= rCutMaxSqr_);
}


bool Foam::interactionLists::testEdgeEdgeDistance
(
    const edge& eI,
    const edge& eJ
) const
{
    const vector& eJs(mesh_.points()[eJ.start()]);
    const vector& eJe(mesh_.points()[eJ.end()]);

    return testEdgeEdgeDistance(eI, eJs, eJe);
}


bool Foam::interactionLists::testPointFaceDistance
(
    const label p,
    const label faceNo
) const
{
    const vector& pointPosition(mesh_.points()[p]);

    return testPointFaceDistance(pointPosition, faceNo);
}


bool Foam::interactionLists::testPointFaceDistance
(
    const label p,
    const referredCell& refCell
) const
{
    const vector& pointPosition(mesh_.points()[p]);

    forAll (refCell.faces(), rCF)
    {
        if
        (
            testPointFaceDistance
            (
                pointPosition,
                refCell.faces()[rCF],
                refCell.vertexPositions(),
                refCell.faceCentres()[rCF],
                refCell.faceAreas()[rCF]
            )
        )
        {
            return true;
        }
    }

    return false;
}


bool Foam::interactionLists::testPointFaceDistance
(
    const vectorList& pointsToTest,
    const label faceNo
) const
{
    forAll(pointsToTest, pTT)
    {
        const vector& p(pointsToTest[pTT]);

        // if any point in the list is in range of the face
        // then the rest do not need to be tested and
        // true can be returned

        if (testPointFaceDistance(p, faceNo))
        {
            return true;
        }
    }

    return false;
}


bool Foam::interactionLists::testPointFaceDistance
(
    const vector& p,
    const label faceNo
) const
{
    const face& faceToTest(mesh_.faces()[faceNo]);

    const vector& faceC(mesh_.faceCentres()[faceNo]);

    const vector& faceA(mesh_.faceAreas()[faceNo]);

    const vectorList& points(mesh_.points());

    return testPointFaceDistance
    (
        p,
        faceToTest,
        points,
        faceC,
        faceA
    );
}


bool Foam::interactionLists::testPointFaceDistance
(
    const vector& p,
    const labelList& faceToTest,
    const vectorList& points,
    const vector& faceC,
    const vector& faceA
) const
{
    vector faceN(faceA/mag(faceA));

    scalar perpDist((p - faceC) & faceN);

    if (magSqr(perpDist) > rCutMaxSqr_)
    {
        return false;
    }

    vector pointOnPlane = (p - faceN * perpDist);

    if (magSqr(faceC - pointOnPlane) < rCutMaxSqr_*1e-8)
    {
        // If pointOnPlane is very close to the face centre
        // then defining the local axes will be inaccurate
        // and it is very likely that pointOnPlane will be
        // inside the face, so return true if the points
        // are in range to be safe

        return (magSqr(pointOnPlane - p) <= rCutMaxSqr_);
    }

    vector xAxis = (faceC - pointOnPlane)/mag(faceC - pointOnPlane);

    vector yAxis =
        ((faceC - pointOnPlane) ^ faceN)
       /mag((faceC - pointOnPlane) ^ faceN);

    List<vector2D> local2DVertices(faceToTest.size());

    forAll(faceToTest, fTT)
    {
        const vector& V(points[faceToTest[fTT]]);

        if (magSqr(V-p) <= rCutMaxSqr_)
        {
            return true;
        }

        local2DVertices[fTT] = vector2D
        (
            ((V - pointOnPlane) & xAxis),
            ((V - pointOnPlane) & yAxis)
        );
    }

    scalar localFaceCx((faceC - pointOnPlane) & xAxis);

    scalar la_valid = -1;

    forAll(local2DVertices, fV)
    {
        const vector2D& va(local2DVertices[fV]);

        const vector2D& vb
        (
            local2DVertices[(fV + 1) % local2DVertices.size()]
        );

        if (mag(vb.y()-va.y()) > SMALL)
        {
            scalar la =
                (
                    va.x() - va.y()*((vb.x() - va.x())/(vb.y() - va.y()))
                )
               /localFaceCx;

            scalar lv = -va.y()/(vb.y() - va.y());


            if (la >= 0 && la <= 1 && lv >= 0 && lv <= 1)
            {
                la_valid = la;

                break;
            }
        }
    }

    if (la_valid < 0)
    {
        // perpendicular point inside face, nearest point is pointOnPlane;
        return (magSqr(pointOnPlane-p) <= rCutMaxSqr_);
    }
    else
    {
        // perpendicular point outside face, nearest point is
        // on edge that generated la_valid;
        return
        (
            magSqr(pointOnPlane + la_valid*(faceC - pointOnPlane) - p)
         <= rCutMaxSqr_
        );
    }

    // if the algorithm hasn't returned anything by now then something has
    // gone wrong.

    FatalErrorIn("interactionLists.C") << nl
        << "point " << p << " to face " << faceToTest
        << " comparison did not find a nearest point"
        << " to be inside or outside face."
        << abort(FatalError);

    return false;
}


bool Foam::interactionLists::testEdgeEdgeDistance
(
    const edge& eI,
    const vector& eJs,
    const vector& eJe
) const
{
    vector a(eI.vec(mesh_.points()));
    vector b(eJe - eJs);

    const vector& eIs(mesh_.points()[eI.start()]);

    vector c(eJs - eIs);

    vector crossab = a ^ b;
    scalar magCrossSqr = magSqr(crossab);

    if (magCrossSqr > VSMALL)
    {
        // If the edges are parallel then a point-face
        // search will pick them up

        scalar s = ((c ^ b) & crossab)/magCrossSqr;
        scalar t = ((c ^ a) & crossab)/magCrossSqr;

        // Check for end points outside of range 0..1
        // If the closest point is outside this range
        // a point-face search will have found it.

        return
        (
            s >= 0
         && s <= 1
         && t >= 0
         && t <= 1
         && magSqr(eIs + a*s - eJs - b*t) <= rCutMaxSqr_
        );
    }

    return false;
}


const Foam::labelList Foam::interactionLists::realCellsInRangeOfSegment
(
    const labelList& segmentFaces,
    const labelList& segmentEdges,
    const labelList& segmentPoints
) const
{
    DynamicList<label> realCellsFoundInRange;

    forAll(segmentFaces, sF)
    {
        const label f = segmentFaces[sF];

        forAll (mesh_.points(), p)
        {
            if (testPointFaceDistance(p, f))
            {
                const labelList& pCells(mesh_.pointCells()[p]);

                forAll(pCells, pC)
                {
                    const label cellI(pCells[pC]);

                    if (findIndex(realCellsFoundInRange, cellI) == -1)
                    {
                        realCellsFoundInRange.append(cellI);
                    }
                }
            }
        }
    }

    forAll(segmentPoints, sP)
    {
        const label p = segmentPoints[sP];

        forAll(mesh_.faces(), f)
        {
            if (testPointFaceDistance(p, f))
            {
                const label cellO(mesh_.faceOwner()[f]);

                if (findIndex(realCellsFoundInRange, cellO) == -1)
                {
                    realCellsFoundInRange.append(cellO);
                }

                if (mesh_.isInternalFace(f))
                {
                    // boundary faces will not have neighbour information

                    const label cellN(mesh_.faceNeighbour()[f]);

                    if (findIndex(realCellsFoundInRange, cellN) == -1)
                    {
                        realCellsFoundInRange.append(cellN);
                    }
                }
            }
        }
    }

    forAll(segmentEdges, sE)
    {
        const edge& eJ(mesh_.edges()[segmentEdges[sE]]);

        forAll (mesh_.edges(), edgeIIndex)
        {
            const edge& eI(mesh_.edges()[edgeIIndex]);

            if (testEdgeEdgeDistance(eI, eJ))
            {
                const labelList& eICells(mesh_.edgeCells()[edgeIIndex]);

                forAll(eICells, eIC)
                {
                    const label cellI(eICells[eIC]);

                    if (findIndex(realCellsFoundInRange, cellI) == -1)
                    {
                        realCellsFoundInRange.append(cellI);
                    }
                }
            }
        }
    }

    return realCellsFoundInRange.shrink();
}


const Foam::labelList Foam::interactionLists::referredCellsInRangeOfSegment
(
    const List<referredCell>& referredInteractionList,
    const labelList& segmentFaces,
    const labelList& segmentEdges,
    const labelList& segmentPoints
) const
{
    DynamicList<label> referredCellsFoundInRange;

    forAll(segmentFaces, sF)
    {
        const label f = segmentFaces[sF];

        forAll(referredInteractionList, rIL)
        {
            const vectorList& refCellPoints
                = referredInteractionList[rIL].vertexPositions();

            if (testPointFaceDistance(refCellPoints, f))
            {
                if (findIndex(referredCellsFoundInRange, rIL) == -1)
                {
                    referredCellsFoundInRange.append(rIL);
                }
            }
        }
    }

    forAll(segmentPoints, sP)
    {
        const label p = segmentPoints[sP];

        forAll(referredInteractionList, rIL)
        {
            const referredCell& refCell(referredInteractionList[rIL]);

            if (testPointFaceDistance(p, refCell))
            {
                if (findIndex(referredCellsFoundInRange, rIL) == -1)
                {
                    referredCellsFoundInRange.append(rIL);
                }
            }
        }
    }

    forAll(segmentEdges, sE)
    {
        const edge& eI(mesh_.edges()[segmentEdges[sE]]);

        forAll(referredInteractionList, rIL)
        {
            const vectorList& refCellPoints
                = referredInteractionList[rIL].vertexPositions();

            const edgeList& refCellEdges
                = referredInteractionList[rIL].edges();

            forAll(refCellEdges, rCE)
            {
                const edge& eJ(refCellEdges[rCE]);

                if
                (
                    testEdgeEdgeDistance
                    (
                        eI,
                        refCellPoints[eJ.start()],
                        refCellPoints[eJ.end()]
                    )
                )
                {
                    if(findIndex(referredCellsFoundInRange, rIL) == -1)
                    {
                        referredCellsFoundInRange.append(rIL);
                    }
                }
            }
        }
    }

    return referredCellsFoundInRange.shrink();
}


// ************************************************************************* //
