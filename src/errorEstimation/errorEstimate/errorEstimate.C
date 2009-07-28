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

#include "errorEstimate.H"
#include "zeroGradientFvPatchField.H"
#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::errorEstimate<Type>::errorBCTypes() const
{
    // Make the boundary condition type list
    // Default types get over-ridden anyway
    wordList ebct
    (
        psi_.boundaryField().size(),
        zeroGradientFvPatchField<Type>::typeName
    );

    forAll (psi_.boundaryField(), patchI)
    {
        if (psi_.boundaryField()[patchI].fixesValue())
        {
            ebct[patchI] = fixedValueFvPatchField<Type>::typeName;
        }
    }

    return ebct;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type>
Foam::errorEstimate<Type>::errorEstimate
(
    const GeometricField<Type, fvPatchField, volMesh>& psi,
    const dimensionSet& ds,
    const Field<Type>& res,
    const scalarField& norm
)
:
    psi_(psi),
    dimensions_(ds),
    residual_(res),
    normFactor_(norm)
{}


// Construct as copy
template<class Type>
Foam::errorEstimate<Type>::errorEstimate(const Foam::errorEstimate<Type>& ee)
:
    refCount(),
    psi_(ee.psi_),
    dimensions_(ee.dimensions_),
    residual_(ee.residual_),
    normFactor_(ee.normFactor_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::errorEstimate<Type>::~errorEstimate()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::errorEstimate<Type>::residual() const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tres
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "residual" + psi_.name(),
                psi_.mesh().time().timeName(),
                psi_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            psi_.dimensions()/dimTime,
            errorBCTypes()
        )
    );

    GeometricField<Type, fvPatchField, volMesh>& res = tres();

    res.internalField() = residual_;
    res.boundaryField() == pTraits<Type>::zero;

    res.correctBoundaryConditions();

    return tres;
}


template<class Type>
Foam::tmp<Foam::volScalarField> Foam::errorEstimate<Type>::normFactor() const
{
    tmp<volScalarField> tnormFactor
    (
        new volScalarField
        (
            IOobject
            (
                "normFactor" + psi_.name(),
                psi_.mesh().time().timeName(),
                psi_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimless/dimTime,
            errorBCTypes()
        )
    );

    volScalarField& normFactor = tnormFactor();

    normFactor.internalField() = normFactor_;
    normFactor.boundaryField() == pTraits<Type>::zero;

    normFactor.correctBoundaryConditions();

    return tnormFactor;
}

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::errorEstimate<Type>::error() const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tresError
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "resError" + psi_.name(),
                psi_.mesh().time().timeName(),
                psi_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            psi_.dimensions(),
            errorBCTypes()
        )
    );

    GeometricField<Type, fvPatchField, volMesh>& resError = tresError();

    resError.internalField() = residual_/normFactor_;
    resError.boundaryField() == pTraits<Type>::zero;

    resError.correctBoundaryConditions();

    return tresError;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::errorEstimate<Type>::operator=(const Foam::errorEstimate<Type>& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "errorEstimate<Type>::operator=(const Foam::errorEstimate<Type>&)"
        )   << "Attempted assignment to self"
            << abort(FatalError);
    }

    if (&psi_ != &(rhs.psi_))
    {
        FatalErrorIn
        (
            "errorEstimate<Type>::operator=(const errorEstimate<Type>&)"
        )   << "different fields"
            << abort(FatalError);
    }

    residual_ = rhs.residual_;
    normFactor_ = rhs.normFactor_;
}


template<class Type>
void Foam::errorEstimate<Type>::operator=(const tmp<errorEstimate<Type> >& teev)
{
    operator=(teev());
    teev.clear();
}


template<class Type>
void Foam::errorEstimate<Type>::negate()
{
    residual_.negate();
}


template<class Type>
void Foam::errorEstimate<Type>::operator+=(const errorEstimate<Type>& eev)
{
    checkMethod(*this, eev, "+=");

    dimensions_ += eev.dimensions_;

    residual_ += eev.residual_;
    normFactor_ += eev.normFactor_;
}


template<class Type>
void Foam::errorEstimate<Type>::operator+=
(
    const tmp<errorEstimate<Type> >& teev
)
{
    operator+=(teev());
    teev.clear();
}


template<class Type>
void Foam::errorEstimate<Type>::operator-=(const errorEstimate<Type>& eev)
{
    checkMethod(*this, eev, "+=");

    dimensions_ -= eev.dimensions_;
    residual_ -= eev.residual_;
    normFactor_ += eev.normFactor_;
}


template<class Type>
void Foam::errorEstimate<Type>::operator-=(const tmp<errorEstimate<Type> >& teev)
{
    operator-=(teev());
    teev.clear();
}


template<class Type>
void Foam::errorEstimate<Type>::operator+=
(
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(*this, su, "+=");
    residual_ -= su.internalField();
}


template<class Type>
void Foam::errorEstimate<Type>::operator+=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    operator+=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::errorEstimate<Type>::operator-=
(
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(*this, su, "-=");
    residual_ += su.internalField();
}


template<class Type>
void Foam::errorEstimate<Type>::operator-=
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    operator-=(tsu());
    tsu.clear();
}


template<class Type>
void Foam::errorEstimate<Type>::operator+=
(
    const dimensioned<Type>& su
)
{
    residual_ -= su;
}


template<class Type>
void Foam::errorEstimate<Type>::operator-=
(
    const dimensioned<Type>& su
)
{
    residual_ += su;
}


template<class Type>
void Foam::errorEstimate<Type>::operator*=
(
    const volScalarField& vsf
)
{
    dimensions_ *= vsf.dimensions();
    residual_ *= vsf.internalField();
    normFactor_ *= vsf.internalField();
}


template<class Type>
void Foam::errorEstimate<Type>::operator*=
(
    const tmp<volScalarField>& tvsf
)
{
    operator*=(tvsf());
    tvsf.clear();
}


template<class Type>
void Foam::errorEstimate<Type>::operator*=
(
    const dimensioned<scalar>& ds
)
{
    dimensions_ *= ds.dimensions();
    residual_ *= ds.value();
    normFactor_ *= ds.value();
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::checkMethod
(
    const errorEstimate<Type>& ee1,
    const errorEstimate<Type>& ee2,
    const char* op
)
{
    if (&ee1.psi() != &ee2.psi())
    {
        FatalErrorIn
        (
            "checkMethod(const errorEstimate<Type>&, "
            "const errorEstimate<Type>&)"
        )   << "incompatible fields for operation "
            << endl << "    "
            << "[" << ee1.psi().name() << "] "
            << op
            << " [" << ee2.psi().name() << "]"
            << abort(FatalError);
    }

    if (dimensionSet::debug && ee1.dimensions() != ee2.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const errorEstimate<Type>&, "
            "const errorEstimate<Type>&)"
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << ee1.psi().name() << ee1.dimensions()/dimVolume << " ] "
            << op
            << " [" << ee2.psi().name() << ee2.dimensions()/dimVolume << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const errorEstimate<Type>& ee,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const char* op
)
{
    if (dimensionSet::debug && ee.dimensions()/dimVolume != vf.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const errorEstimate<Type>&, "
            "const GeometricField<Type, fvPatchField, volMesh>&)"
        )   <<  "incompatible dimensions for operation "
            << endl << "    "
            << "[" << ee.psi().name() << ee.dimensions()/dimVolume << " ] "
            << op
            << " [" << vf.name() << vf.dimensions() << " ]"
            << abort(FatalError);
    }
}


template<class Type>
void Foam::checkMethod
(
    const errorEstimate<Type>& ee,
    const dimensioned<Type>& dt,
    const char* op
)
{
    if (dimensionSet::debug && ee.dimensions()/dimVolume != dt.dimensions())
    {
        FatalErrorIn
        (
            "checkMethod(const errorEstimate<Type>&, const dimensioned<Type>&)"
        )   << "incompatible dimensions for operation "
            << endl << "    "
            << "[" << ee.psi().name() << ee.dimensions()/dimVolume << " ] "
            << op
            << " [" << dt.name() << dt.dimensions() << " ]"
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const errorEstimate<Type>& A,
    const errorEstimate<Type>& B
)
{
    checkMethod(A, B, "+");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC() += B;
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<errorEstimate<Type> >& tA,
    const errorEstimate<Type>& B
)
{
    checkMethod(tA(), B, "+");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC() += B;
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const errorEstimate<Type>& A,
    const tmp<errorEstimate<Type> >& tB
)
{
    checkMethod(A, tB(), "+");
    tmp<errorEstimate<Type> > tC(tB.ptr());
    tC() += A;
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<errorEstimate<Type> >& tA,
    const tmp<errorEstimate<Type> >& tB
)
{
    checkMethod(tA(), tB(), "+");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC() += tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const errorEstimate<Type>& A
)
{
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().negate();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<errorEstimate<Type> >& tA
)
{
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().negate();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const errorEstimate<Type>& A,
    const errorEstimate<Type>& B
)
{
    checkMethod(A, B, "-");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC() -= B;
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<errorEstimate<Type> >& tA,
    const errorEstimate<Type>& B
)
{
    checkMethod(tA(), B, "-");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC() -= B;
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const errorEstimate<Type>& A,
    const tmp<errorEstimate<Type> >& tB
)
{
    checkMethod(A, tB(), "-");
    tmp<errorEstimate<Type> > tC(tB.ptr());
    tC() -= A;
    tC().negate();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<errorEstimate<Type> >& tA,
    const tmp<errorEstimate<Type> >& tB
)
{
    checkMethod(tA(), tB(), "-");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC() -= tB();
    tB.clear();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const errorEstimate<Type>& A,
    const errorEstimate<Type>& B
)
{
    checkMethod(A, B, "==");
    return (A - B);
}


template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const tmp<errorEstimate<Type> >& tA,
    const errorEstimate<Type>& B
)
{
    checkMethod(tA(), B, "==");
    return (tA - B);
}


template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const errorEstimate<Type>& A,
    const tmp<errorEstimate<Type> >& tB
)
{
    checkMethod(A, tB(), "==");
    return (A - tB);
}


template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const tmp<errorEstimate<Type> >& tA,
    const tmp<errorEstimate<Type> >& tB
)
{
    checkMethod(tA(), tB(), "==");
    return (tA - tB);
}


template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const errorEstimate<Type>& A,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(A, su, "+");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().res() -= su.internalField();
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<errorEstimate<Type> >& tA,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().res() -= su.internalField();
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const errorEstimate<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "+");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().res() -= tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<errorEstimate<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().res() -= tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const GeometricField<Type, fvPatchField, volMesh>& su,
    const errorEstimate<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().res() -= su.internalField();
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const GeometricField<Type, fvPatchField, volMesh>& su,
    const tmp<errorEstimate<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().res() -= su.internalField();
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const errorEstimate<Type>& A
)
{
    checkMethod(A, tsu(), "+");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().res() -= tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const tmp<errorEstimate<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "+");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().res() -= tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const errorEstimate<Type>& A,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(A, su, "-");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().res() += su.internalField();
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<errorEstimate<Type> >& tA,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().res() += su.internalField();
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const errorEstimate<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "-");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().res() += tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<errorEstimate<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().res() += tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const GeometricField<Type, fvPatchField, volMesh>& su,
    const errorEstimate<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().negate();
    tC().res() -= su.internalField();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const GeometricField<Type, fvPatchField, volMesh>& su,
    const tmp<errorEstimate<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().negate();
    tC().res() -= su.internalField();
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const errorEstimate<Type>& A
)
{
    checkMethod(A, tsu(), "-");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().negate();
    tC().res() -= tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu,
    const tmp<errorEstimate<Type> >& tA
)
{
    checkMethod(tA(), tsu(), "-");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().negate();
    tC().res() -= tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const errorEstimate<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "+");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().res() -= su.value();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const tmp<errorEstimate<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "+");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().res() -= su.value();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const dimensioned<Type>& su,
    const errorEstimate<Type>& A
)
{
    checkMethod(A, su, "+");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().res() -= su.value();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator+
(
    const dimensioned<Type>& su,
    const tmp<errorEstimate<Type> >& tA
)
{
    checkMethod(tA(), su, "+");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().res() -= su.value();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const errorEstimate<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "-");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().res() += su.value();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const tmp<errorEstimate<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "-");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().res() += su.value();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const dimensioned<Type>& su,
    const errorEstimate<Type>& A
)
{
    checkMethod(A, su, "-");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().negate();
    tC().res() -= su.value();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator-
(
    const dimensioned<Type>& su,
    const tmp<errorEstimate<Type> >& tA
)
{
    checkMethod(tA(), su, "-");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().negate();
    tC().res() -= su.value();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const errorEstimate<Type>& A,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(A, su, "==");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().res() += su.internalField();
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const tmp<errorEstimate<Type> >& tA,
    const GeometricField<Type, fvPatchField, volMesh>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().res() += su.internalField();
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const errorEstimate<Type>& A,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(A, tsu(), "==");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().res() += tsu().internalField();
    tsu.clear();
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const tmp<errorEstimate<Type> >& tA,
    const tmp<GeometricField<Type, fvPatchField, volMesh> >& tsu
)
{
    checkMethod(tA(), tsu(), "==");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().res() += tsu().internalField();
    tsu.clear();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const errorEstimate<Type>& A,
    const dimensioned<Type>& su
)
{
    checkMethod(A, su, "==");
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC().res() += su.value();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator==
(
    const tmp<errorEstimate<Type> >& tA,
    const dimensioned<Type>& su
)
{
    checkMethod(tA(), su, "==");
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC().res() += su.value();
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator*
(
    const volScalarField& vsf,
    const errorEstimate<Type>& A
)
{
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC() *= vsf;
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator*
(
    const tmp<volScalarField>& tvsf,
    const errorEstimate<Type>& A
)
{
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC() *= tvsf;
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator*
(
    const volScalarField& vsf,
    const tmp<errorEstimate<Type> >& tA
)
{
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC() *= vsf;
    return tC;
}

template<class Type>
tmp<errorEstimate<Type> > operator*
(
    const tmp<volScalarField>& tvsf,
    const tmp<errorEstimate<Type> >& tA
)
{
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC() *= tvsf;
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator*
(
    const dimensioned<scalar>& ds,
    const errorEstimate<Type>& A
)
{
    tmp<errorEstimate<Type> > tC(new errorEstimate<Type>(A));
    tC() *= ds;
    return tC;
}


template<class Type>
tmp<errorEstimate<Type> > operator*
(
    const dimensioned<scalar>& ds,
    const tmp<errorEstimate<Type> >& tA
)
{
    tmp<errorEstimate<Type> > tC(tA.ptr());
    tC() *= ds;
    return tC;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
