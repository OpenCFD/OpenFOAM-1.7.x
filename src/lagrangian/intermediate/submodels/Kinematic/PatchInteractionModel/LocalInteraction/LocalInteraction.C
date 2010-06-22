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

#include "LocalInteraction.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template <class CloudType>
Foam::label Foam::LocalInteraction<CloudType>::applyToPatch
(
    const label globalPatchI
) const
{
    forAll(patchIds_, patchI)
    {
        if (patchIds_[patchI] == globalPatchI)
        {
            return patchI;
        }
    }

    return -1;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::LocalInteraction<CloudType>::LocalInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    patchData_(this->coeffDict().lookup("patches")),
    patchIds_(patchData_.size())
{
    const polyMesh& mesh = cloud.mesh();
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    // check that user patches are valid region patches
    forAll(patchData_, patchI)
    {
        const word& patchName = patchData_[patchI].patchName();
        patchIds_[patchI] = bMesh.findPatchID(patchName);
        if (patchIds_[patchI] < 0)
        {
            FatalErrorIn("LocalInteraction(const dictionary&, CloudType&)")
                << "Patch " << patchName << " not found. Available patches "
                << "are: " << bMesh.names() << nl << exit(FatalError);
        }
    }

    // check that all walls are specified
    DynamicList<word> badWalls;
    forAll(bMesh, patchI)
    {
        if
        (
            isA<wallPolyPatch>(bMesh[patchI])
         && applyToPatch(bMesh[patchI].index()) < 0
        )
        {
            badWalls.append(bMesh[patchI].name());
        }
    }

    if (badWalls.size() > 0)
    {
        FatalErrorIn("LocalInteraction(const dictionary&, CloudType&)")
            << "All wall patches must be specified when employing local patch "
            << "interaction. Please specify data for patches:" << nl
            << badWalls << nl << exit(FatalError);
    }

    // check that interactions are valid/specified
    forAll(patchData_, patchI)
    {
        const word& interactionTypeName =
            patchData_[patchI].interactionTypeName();
        const typename PatchInteractionModel<CloudType>::interactionType& it =
            this->wordToInteractionType(interactionTypeName);

        if (it == PatchInteractionModel<CloudType>::itOther)
        {
            const word& patchName = patchData_[patchI].patchName();
            FatalErrorIn("LocalInteraction(const dictionary&, CloudType&)")
                << "Unknown patch interaction type "
                << interactionTypeName << " for patch " << patchName
                << ". Valid selections are:"
                << this->PatchInteractionModel<CloudType>::interactionTypeNames_
                << nl << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::LocalInteraction<CloudType>::~LocalInteraction()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::LocalInteraction<CloudType>::active() const
{
    return true;
}


template <class CloudType>
bool Foam::LocalInteraction<CloudType>::correct
(
    const polyPatch& pp,
    const label faceId,
    bool& keepParticle,
    bool& active,
    vector& U
) const
{
    label patchI = applyToPatch(pp.index());

    if (patchI >= 0)
    {
        typename PatchInteractionModel<CloudType>::interactionType it =
            this->wordToInteractionType
            (
                patchData_[patchI].interactionTypeName()
            );

        switch (it)
        {
            case PatchInteractionModel<CloudType>::itEscape:
            {
                keepParticle = false;
                active = false;
                U = vector::zero;
                break;
            }
            case PatchInteractionModel<CloudType>::itStick:
            {
                keepParticle = true;
                active = false;
                U = vector::zero;
                break;
            }
            case PatchInteractionModel<CloudType>::itRebound:
            {
                keepParticle = true;
                active = true;

                vector nw = pp.faceAreas()[pp.whichFace(faceId)];
                nw /= mag(nw);

                scalar Un = U & nw;
                vector Ut = U - Un*nw;

                if (Un > 0)
                {
                    U -= (1.0 + patchData_[patchI].e())*Un*nw;
                }

                U -= patchData_[patchI].mu()*Ut;

                break;
            }
            default:
            {
                FatalErrorIn
                (
                    "bool LocalInteraction<CloudType>::correct"
                    "("
                        "const polyPatch&, "
                        "const label, "
                        "bool&, "
                        "vector&"
                    ") const"
                )   << "Unknown interaction type "
                    << patchData_[patchI].interactionTypeName()
                    << "(" << it << ") for patch "
                    << patchData_[patchI].patchName()
                    << ". Valid selections are:" << this->interactionTypeNames_
                    << endl << abort(FatalError);
            }
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
