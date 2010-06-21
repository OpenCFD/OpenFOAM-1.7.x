/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "fv.H"
#include "snGradScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

template<class Type>
tmp<snGradScheme<Type> > snGradScheme<Type>::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (fv::debug)
    {
        Info<< "snGradScheme<Type>::New(const fvMesh&, Istream&)"
               " : constructing snGradScheme<Type>"
            << endl;
    }

    if (schemeData.eof())
    {
        FatalIOErrorIn
        (
            "snGradScheme<Type>::New(const fvMesh&, Istream&)",
            schemeData
        )   << "Discretisation scheme not specified"
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    word schemeName(schemeData);

    typename MeshConstructorTable::iterator constructorIter =
        MeshConstructorTablePtr_->find(schemeName);

    if (constructorIter == MeshConstructorTablePtr_->end())
    {
        FatalIOErrorIn
        (
            "snGradScheme<Type>::New(const fvMesh&, Istream&)",
            schemeData
        )   << "Unknown discretisation scheme " << schemeName
            << endl << endl
            << "Valid schemes are :" << endl
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return constructorIter()(mesh, schemeData);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
snGradScheme<Type>::~snGradScheme()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
snGradScheme<Type>::snGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const tmp<surfaceScalarField>& tdeltaCoeffs,
    const word& snGradName
)
{
    const fvMesh& mesh = vf.mesh();

    // construct GeometricField<Type, fvsPatchField, surfaceMesh>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tssf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                snGradName + "("+vf.name()+')',
                vf.instance(),
                vf.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vf.dimensions()*tdeltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& ssf = tssf();

    // set reference to difference factors array
    const scalarField& deltaCoeffs = tdeltaCoeffs().internalField();

    // owner/neighbour addressing
    const unallocLabelList& owner = mesh.owner();
    const unallocLabelList& neighbour = mesh.neighbour();

    forAll(owner, faceI)
    {
        ssf[faceI] =
            deltaCoeffs[faceI]*(vf[neighbour[faceI]] - vf[owner[faceI]]);
    }

    forAll(vf.boundaryField(), patchI)
    {
        ssf.boundaryField()[patchI] = vf.boundaryField()[patchI].snGrad();
    }

    return tssf;
}


//- Return the face-snGrad of the given cell field
//  with explicit correction
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
snGradScheme<Type>::snGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
        = snGrad(vf, deltaCoeffs(vf));

    if (corrected())
    {
        tsf() += correction(vf);
    }

    return tsf;
}


//- Return the face-snGrad of the given cell field
//  with explicit correction
template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
snGradScheme<Type>::snGrad
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf
) const
{
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tinterpVf
        = snGrad(tvf());
    tvf.clear();
    return tinterpVf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
