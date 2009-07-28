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

#include "meshToMesh.H"
#include "volFields.H"
#include "interpolationCellPoint.H"
#include "SubField.H"
#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void meshToMesh::mapField
(
    Field<Type>& toF,
    const Field<Type>& fromVf,
    const labelList& adr
) const
{
    // Direct mapping of nearest-cell values

    forAll(toF, celli)
    {
        if (adr[celli] != -1)
        {
            toF[celli] = fromVf[adr[celli]];
        }
    }

    //toF.map(fromVf, adr);
}


template<class Type>
void meshToMesh::interpolateField
(
    Field<Type>& toF,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf,
    const labelList& adr,
    const scalarListList& weights
) const
{
    // Inverse distance weighted interpolation

    // get reference to cellCells
    const labelListList& cc = fromMesh_.cellCells();

    forAll (toF, celli)
    {
        if (adr[celli] != -1)
        {
            const labelList& neighbours = cc[adr[celli]];
            const scalarList& w = weights[celli];
            
            toF[celli] = fromVf[adr[celli]]*w[0];

            for (label ni = 1; ni < w.size(); ni++)
            {
                toF[celli] += fromVf[neighbours[ni - 1]]*w[ni];
            }
        }
    }
}


template<class Type>
void meshToMesh::interpolateField
(
    Field<Type>& toF,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf,
    const labelList& adr,
    const vectorField& centres
) const
{
    // Cell-Point interpolation
    interpolationCellPoint<Type> interpolator(fromVf);

    forAll (toF, celli)
    {
        if (adr[celli] != -1)
        {
            toF[celli] = interpolator.interpolate
            (
                centres[celli],
                adr[celli]
            );
        }
    }
}


template<class Type>
void meshToMesh::interpolateInternalField
(
    Field<Type>& toF,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf,
    meshToMesh::order ord
) const
{
    if (fromVf.mesh() != fromMesh_)
    {
        FatalErrorIn
        (
            "meshToMesh::interpolateInternalField(Field<Type>& toF, "
            "const GeometricField<Type, fvPatchField, volMesh>& fromVf, "
            "meshToMesh::order ord) const"
        )   << "the argument field does not correspond to the right mesh. "
            << "Field size: " << fromVf.size()
            << " mesh size: " << fromMesh_.nCells()
            << exit(FatalError);
    }

    if (toF.size() != toMesh_.nCells())
    {
        FatalErrorIn
        (
            "meshToMesh::interpolateInternalField(Field<Type>& toF, "
            "const GeometricField<Type, fvPatchField, volMesh>& fromVf, "
            "meshToMesh::order ord) const"
        )   << "the argument field does not correspond to the right mesh. "
            << "Field size: " << toF.size()
            << " mesh size: " << toMesh_.nCells()
            << exit(FatalError);
    }

    switch(ord)
    {
        case MAP:
            mapField(toF, fromVf, cellAddressing_);
        break;

        case INTERPOLATE:
            interpolateField
            (
                toF, 
                fromVf,
                cellAddressing_,
                inverseDistanceWeights()
            );
        break;

        case CELL_POINT_INTERPOLATE:
            interpolateField
            (
                toF, 
                fromVf,
                cellAddressing_,
                toMesh_.cellCentres()
            );
        break;

        default:
            FatalErrorIn
            (
                "meshToMesh::interpolateInternalField(Field<Type>& toF, "
                "const GeometricField<Type, fvPatchField, volMesh>& fromVf, "
                "meshToMesh::order ord) const"
            )   << "unknown interpolation scheme " << ord
                << exit(FatalError);
    }
}


template<class Type>
void meshToMesh::interpolateInternalField
(
    Field<Type>& toF,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfromVf,
    meshToMesh::order ord
) const
{
    interpolateInternalField(toF, tfromVf(), ord);
    tfromVf.clear();
}


template<class Type>
void meshToMesh::interpolate
(
    GeometricField<Type, fvPatchField, volMesh>& toVf,
    const GeometricField<Type, fvPatchField, volMesh>& fromVf,
    meshToMesh::order ord
) const
{
    interpolateInternalField(toVf, fromVf, ord);

    forAll (toMesh_.boundaryMesh(), patchi)
    {
        const fvPatch& toPatch = toMesh_.boundary()[patchi];

        if (cuttingPatches_.found(toPatch.name()))
        {
            switch(ord)
            {
                case MAP:
                    mapField
                    (
                        toVf.boundaryField()[patchi],
                        fromVf,
                        boundaryAddressing_[patchi]
                    );
                break;

                case INTERPOLATE:
                    interpolateField
                    (
                        toVf.boundaryField()[patchi],
                        fromVf,
                        boundaryAddressing_[patchi],
                        toPatch.Cf()
                    );
                break;

                case CELL_POINT_INTERPOLATE:
                    interpolateField
                    (
                        toVf.boundaryField()[patchi],
                        fromVf,
                        boundaryAddressing_[patchi],
                        toPatch.Cf()
                    );
                break;

                default:
                    FatalErrorIn
                    (
                        "meshToMesh::interpolate("
                        "GeometricField<Type, fvPatchField, volMesh>& toVf, "
                        "const GeometricField<Type, fvPatchField, volMesh>& "
                        "fromVf, meshToMesh::order ord) const"
                    )   << "unknown interpolation scheme " << ord
                        << exit(FatalError);
            }

            if (isA<mixedFvPatchField<Type> >(toVf.boundaryField()[patchi]))
            {
                refCast<mixedFvPatchField<Type> >
                (
                    toVf.boundaryField()[patchi]
                ).refValue() = toVf.boundaryField()[patchi];
            }
        }
        else if 
        (
            patchMap_.found(toPatch.name())
         && fromMeshPatches_.found(patchMap_.find(toPatch.name())())
        )
        {
            /*
            toVf.boundaryField()[patchi].map
            (
                fromVf.boundaryField()
                [
                    fromMeshPatches_.find(patchMap_.find(toPatch.name())())()
                ],
                boundaryAddressing_[patchi]
            );
            */

            mapField
            (
                toVf.boundaryField()[patchi],
                fromVf.boundaryField()
                [
                    fromMeshPatches_.find(patchMap_.find(toPatch.name())())()
                ],
                boundaryAddressing_[patchi]
            );
        }
    }
}


template<class Type>
void meshToMesh::interpolate
(
    GeometricField<Type, fvPatchField, volMesh>& toVf,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfromVf,
    meshToMesh::order ord
) const
{
    interpolate(toVf, tfromVf(), ord);
    tfromVf.clear();
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > meshToMesh::interpolate
(
    const GeometricField<Type, fvPatchField, volMesh>& fromVf,
    meshToMesh::order ord
) const
{
    // Create and map the internal-field values
    Field<Type> internalField(toMesh_.nCells());
    interpolateInternalField(internalField, fromVf, ord);

    // check whether both meshes have got the same number
    // of boundary patches
    if (fromMesh_.boundary().size() != toMesh_.boundary().size())
    {
        FatalErrorIn
        (
            "meshToMesh::interpolate"
            "(const GeometricField<Type, fvPatchField, volMesh>& fromVf,"
            "meshToMesh::order ord) const"
        )   << "Incompatible meshes: different number of boundaries, "
               "only internal field may be interpolated"
            << exit(FatalError);
    }

    // Create and map the patch field values
    PtrList<fvPatchField<Type> > patchFields
    (
        boundaryAddressing_.size()
    );

    forAll (boundaryAddressing_, patchI)
    {
        patchFields.set
        (
            patchI,
            fvPatchField<Type>::New
            (
                fromVf.boundaryField()[patchI],
                toMesh_.boundary()[patchI],
                DimensionedField<Type, volMesh>::null(),
                patchFieldInterpolator
                (
                    boundaryAddressing_[patchI]
                )
            )
        );
    }


    // Create the complete field from the pieces
    tmp<GeometricField<Type, fvPatchField, volMesh> > ttoF
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "interpolated(" + fromVf.name() + ')',
                toMesh_.time().timeName(),
                toMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            toMesh_,
            fromVf.dimensions(),
            internalField,
            patchFields
        )
    );

    return ttoF;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > meshToMesh::interpolate
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tfromVf,
    meshToMesh::order ord
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tint =
        interpolate(tfromVf(), ord);
    tfromVf.clear();

    return tint;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
