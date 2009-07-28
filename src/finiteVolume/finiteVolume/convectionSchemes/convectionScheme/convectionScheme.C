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

Description
    Abstract base class for finite volume calculus convection schemes.

\*---------------------------------------------------------------------------*/

#include "fv.H"
#include "HashTable.H"
#include "linear.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
convectionScheme<Type>::convectionScheme(const convectionScheme& cs)
:
    refCount(),
    mesh_(cs.mesh_)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<convectionScheme<Type> > convectionScheme<Type>::New
(
    const fvMesh& mesh,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        Info<< "convectionScheme<Type>::New"
               "(const fvMesh&, const surfaceScalarField&, Istream&) : "
               "constructing convectionScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "convectionScheme<Type>::New"
            "(const fvMesh&, const surfaceScalarField&, Istream&)",
            schemeData
        )   << "Convection scheme not specified" << endl << endl
            << "Valid convection schemes are :" << endl
            << IstreamConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    word schemeName(schemeData);

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(schemeName);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "convectionScheme<Type>::New"
            "(const fvMesh&, const surfaceScalarField&, Istream&)",
            schemeData
        )   << "unknown convection scheme " << schemeName << endl << endl
            << "Valid convection schemes are :" << endl
            << IstreamConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, faceFlux, schemeData);
}


template<class Type>
tmp<convectionScheme<Type> > convectionScheme<Type>::New
(
    const fvMesh& mesh,
    const typename multivariateSurfaceInterpolationScheme<Type>::
        fieldTable& fields,
    const surfaceScalarField& faceFlux,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        Info<< "convectionScheme<Type>::New"
               "(const fvMesh&, "
               "const typename multivariateSurfaceInterpolationScheme<Type>"
               "::fieldTable&, const surfaceScalarField&, Istream&) : "
               "constructing convectionScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "convectionScheme<Type>::New"
               "(const fvMesh&, "
               "const typename multivariateSurfaceInterpolationScheme<Type>"
               "::fieldTable&, const surfaceScalarField&, Istream&)",
            schemeData
        )   << "Convection scheme not specified" << endl << endl
            << "Valid convection schemes are :" << endl
            << MultivariateConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    word schemeName(schemeData);

    typename MultivariateConstructorTable::iterator cstrIter =
        MultivariateConstructorTablePtr_->find(schemeName);

    if (cstrIter == MultivariateConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "convectionScheme<Type>::New"
            "(const fvMesh&, "
            "const typename multivariateSurfaceInterpolationScheme<Type>"
            "::fieldTable&, const surfaceScalarField&, Istream&)",
            schemeData
        )   << "unknown convection scheme " << schemeName << endl << endl
            << "Valid convection schemes are :" << endl
            << MultivariateConstructorTablePtr_->toc()
            << exit(FatalIOError);
    }

    return cstrIter()(mesh, fields, faceFlux, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
convectionScheme<Type>::~convectionScheme()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void convectionScheme<Type>::operator=(const convectionScheme<Type>& cs)
{
    if (this == &cs)
    {
        FatalErrorIn
        (
            "convectionScheme<Type>::operator=(const convectionScheme<Type>&)"
        )   << "attempted assignment to self"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
