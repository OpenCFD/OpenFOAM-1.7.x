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

\*---------------------------------------------------------------------------*/

#include "PatchPostProcessing.H"
#include "IOPtrList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <class CloudType>
Foam::label Foam::PatchPostProcessing<CloudType>::applyToPatch
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


// * * * * * * * * * * * * * protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::PatchPostProcessing<CloudType>::write()
{
    forAll(patchData_, patchI)
    {
        IOPtrList<parcelType> postObject
        (
            IOobject
            (
                patchNames_[patchI] + ".post",
                this->owner().time().timeName(),
                "postProcessing",
                this->owner(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            patchData_[patchI].size()
        );

        forAll(postObject, ptrI)
        {
            postObject.set(ptrI, patchData_[patchI][ptrI].ptr());
        }

        postObject.note() = parcelType::propHeader;

        postObject.writeObject
        (
            IOstream::ASCII,
            IOstream::currentVersion,
            mesh_.time().writeCompression()
        );

        patchData_[patchI].clearStorage();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchPostProcessing<CloudType>::PatchPostProcessing
(
    const dictionary& dict,
    CloudType& owner
)
:
    PostProcessingModel<CloudType>(dict, owner, typeName),
    mesh_(owner.mesh()),
    patchNames_(this->coeffDict().lookup("patches")),
    patchData_(patchNames_.size()),
    patchIds_(patchNames_.size())
{
    forAll(patchNames_, patchI)
    {
        label id = mesh_.boundaryMesh().findPatchID(patchNames_[patchI]);
        if (id < 0)
        {
            FatalErrorIn
            (
                "PatchPostProcessing<CloudType>::PatchPostProcessing"
                "("
                    "const dictionary&, "
                    "CloudType& owner"
                ")"
            )<< "Requested patch " << patchNames_[patchI] << " not found" << nl
             << "Available patches are: " << mesh_.boundaryMesh().names() << nl
             << exit(FatalError);
        }
        patchIds_[patchI] = id;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PatchPostProcessing<CloudType>::~PatchPostProcessing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::PatchPostProcessing<CloudType>::active() const
{
    return true;
}


template<class CloudType>
void Foam::PatchPostProcessing<CloudType>::postPatch
(
    const parcelType& p,
    const label patchI
)
{
    label localPatchI = applyToPatch(patchI);
    if (localPatchI >= 0 && patchData_[localPatchI].size() < maxStoredParcels_)
    {
        patchData_[localPatchI].append(p.clone());
    }
}


// ************************************************************************* //
