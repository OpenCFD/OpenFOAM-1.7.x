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

#include "motionSmoother.H"
#include "meshTools.H"
#include "processorPointPatchFields.H"
#include "globalPointPatchFields.H"
#include "pointConstraint.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
void Foam::motionSmoother::checkConstraints
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> FldType;

    const polyMesh& mesh = pf.mesh();

    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    // first count the total number of patch-patch points

    label nPatchPatchPoints = 0;

    forAll(bm, patchi)
    {
        if(!isA<emptyPolyPatch>(bm[patchi]))
        {
            nPatchPatchPoints += bm[patchi].boundaryPoints().size();
        }
    }


    typename FldType::GeometricBoundaryField& bFld = pf.boundaryField();


    // Evaluate in reverse order

    forAllReverse(bFld, patchi)
    {
        bFld[patchi].initEvaluate(Pstream::blocking);   // buffered
    }

    forAllReverse(bFld, patchi)
    {
        bFld[patchi].evaluate(Pstream::blocking);
    }


    // Save the values

    Field<Type> boundaryPointValues(nPatchPatchPoints);
    nPatchPatchPoints = 0;

    forAll(bm, patchi)
    {
        if(!isA<emptyPolyPatch>(bm[patchi]))
        {
            const labelList& bp = bm[patchi].boundaryPoints();
            const labelList& meshPoints = bm[patchi].meshPoints();

            forAll(bp, pointi)
            {
                label ppp = meshPoints[bp[pointi]];
                boundaryPointValues[nPatchPatchPoints++] = pf[ppp];
            }
        }
    }
    

    // Forward evaluation

    bFld.evaluate();


    // Check

    nPatchPatchPoints = 0;

    forAll(bm, patchi)
    {
        if(!isA<emptyPolyPatch>(bm[patchi]))
        {
            const labelList& bp = bm[patchi].boundaryPoints();
            const labelList& meshPoints = bm[patchi].meshPoints();

            forAll(bp, pointi)
            {
                label ppp = meshPoints[bp[pointi]];

                const Type& savedVal = boundaryPointValues[nPatchPatchPoints++];

                if (savedVal != pf[ppp])
                {
                    FatalErrorIn
                    (
                        "motionSmoother::checkConstraints"
                        "(GeometricField<Type, pointPatchField, pointMesh>&)"
                    )   << "Patch fields are not consistent on mesh point "
                        << ppp << " coordinate " << mesh.points()[ppp]
                        << " at patch " << bm[patchi].name() << '.'
                        << endl
                        << "Reverse evaluation gives value " << savedVal
                        << " , forward evaluation gives value " << pf[ppp]
                        << abort(FatalError);
                }
            }
        }
    }
}


template<class Type>
void Foam::motionSmoother::applyCornerConstraints
(
    GeometricField<Type, pointPatchField, pointMesh>& pf
) const
{
    forAll(patchPatchPointConstraintPoints_, pointi)
    {
        pf[patchPatchPointConstraintPoints_[pointi]] = transform
        (
            patchPatchPointConstraintTensors_[pointi],
            pf[patchPatchPointConstraintPoints_[pointi]]
        );
    }
}


// Average of connected points.
template <class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh> >
 Foam::motionSmoother::avg
(
    const GeometricField<Type, pointPatchField, pointMesh>& fld,
    const scalarField& edgeWeight,
    const bool separation
) const
{
    tmp<GeometricField<Type, pointPatchField, pointMesh> > tres
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                "avg("+fld.name()+')',
                fld.time().timeName(),
                fld.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fld.mesh(),
            dimensioned<Type>("zero", fld.dimensions(), pTraits<Type>::zero)
        )
    );
    GeometricField<Type, pointPatchField, pointMesh>& res = tres();

    const polyMesh& mesh = fld.mesh()();


    // Sum local weighted values and weights
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Note: on coupled edges use only one edge (through isMasterEdge)
    // This is done so coupled edges do not get counted double.

    scalarField sumWeight(mesh.nPoints(), 0.0);

    const edgeList& edges = mesh.edges();

    forAll(edges, edgeI)
    {
        if (isMasterEdge_.get(edgeI) == 1)
        {
            const edge& e = edges[edgeI];
            const scalar w = edgeWeight[edgeI];

            res[e[0]] += w*fld[e[1]];
            sumWeight[e[0]] += w;

            res[e[1]] += w*fld[e[0]];
            sumWeight[e[1]] += w;
        }
    }


    // Add coupled contributions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    syncTools::syncPointList
    (
        mesh,
        res,
        plusEqOp<Type>(),
        pTraits<Type>::zero,    // null value
        separation              // separation
    );

    syncTools::syncPointList
    (
        mesh,
        sumWeight,
        plusEqOp<scalar>(),
        scalar(0),              // null value
        false                   // separation
    );


    // Average
    // ~~~~~~~

    forAll(res, pointI)
    {
        if (mag(sumWeight[pointI]) < VSMALL)
        {
            // Unconnected point. Take over original value
            res[pointI] = fld[pointI];
        }
        else
        {
            res[pointI] /= sumWeight[pointI];
        }        
    }

    res.correctBoundaryConditions();
    applyCornerConstraints(res);

    return tres;
}


// smooth field (point-jacobi)
template <class Type>
void Foam::motionSmoother::smooth
(
    const GeometricField<Type, pointPatchField, pointMesh>& fld,
    const scalarField& edgeWeight,
    const bool separation,
    GeometricField<Type, pointPatchField, pointMesh>& newFld
) const
{
    tmp<pointVectorField> tavgFld = avg
    (
        fld,
        edgeWeight, // weighting
        separation  // whether to apply separation vector
    );
    const pointVectorField& avgFld = tavgFld();

    forAll(fld, pointI)
    {
        if (isInternalPoint(pointI))
        {
            newFld[pointI] = 0.5*fld[pointI] + 0.5*avgFld[pointI];
        }
    }

    newFld.correctBoundaryConditions();
    applyCornerConstraints(newFld);
}


//- Test synchronisation of pointField
template<class Type, class CombineOp>
void Foam::motionSmoother::testSyncField
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const Type& zero,
    const bool separation,
    const scalar maxMag
) const
{
    if (debug)
    {
        Pout<< "testSyncField : testing synchronisation of Field<Type>."
            << endl;
    }

    Field<Type> syncedFld(fld);

    syncTools::syncPointList
    (
        mesh_,
        syncedFld,
        cop,
        zero,       // null value
        separation  // separation
    );

    forAll(syncedFld, i)
    {
        if (mag(syncedFld[i] - fld[i]) > maxMag)
        {
            FatalErrorIn
            (
                "motionSmoother::testSyncField"
                "(const Field<Type>&, const CombineOp&"
                ", const Type&, const bool)"
            )   << "On element " << i << " value:" << fld[i]
                << " synchronised value:" << syncedFld[i]
                << abort(FatalError);
        }
    }
}


// ************************************************************************* //
