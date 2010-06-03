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

#include "Polynomial.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial()
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0.0)
{}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial(const word& name, Istream& is)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(),
    logActive_(false),
    logCoeff_(0.0)
{
    word isName(is);

    if (isName != name)
    {
        FatalErrorIn
        (
            "Polynomial<PolySize>::Polynomial(const word&, Istream&)"
        )   << "Expected polynomial name " << name << " but read " << isName
            << nl << exit(FatalError);
    }

    VectorSpace<Polynomial<PolySize>, scalar, PolySize>::
        operator=(VectorSpace<Polynomial<PolySize>, scalar, PolySize>(is));

    if (this->size() == 0)
    {
        FatalErrorIn
        (
            "Polynomial<PolySize>::Polynomial(const word&, Istream&)"
        )   << "Polynomial coefficients for entry " << isName
            << " are invalid (empty)" << nl << exit(FatalError);
    }
}


template<int PolySize>
Foam::Polynomial<PolySize>::Polynomial
(
    const Polynomial<PolySize>& poly
)
:
    VectorSpace<Polynomial<PolySize>, scalar, PolySize>(poly),
    logActive_(poly.logActive_),
    logCoeff_(poly.logCoeff_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<int PolySize>
bool& Foam::Polynomial<PolySize>::logActive()
{
    return logActive_;
}


template<int PolySize>
Foam::scalar& Foam::Polynomial<PolySize>::logCoeff()
{
    return logCoeff_;
}


template<int PolySize>
Foam::scalar Foam::Polynomial<PolySize>::evaluate(const scalar x) const
{
    scalar y = this->v_[0];

    for (label i=1; i<PolySize; i++)
    {
        y += this->v_[i]*pow(x, i);
    }

    if (logActive_)
    {
        y += logCoeff_*log(x);
    }

    return y;
}


template<int PolySize>
Foam::scalar Foam::Polynomial<PolySize>::integrateLimits
(
    const scalar x1,
    const scalar x2
) const
{
    if (logActive_)
    {
        FatalErrorIn
        (
            "scalar Polynomial<PolySize>::integrateLimits"
            "("
                "const scalar, "
                "const scalar"
            ") const"
        )   << "Cannot integrate polynomial with logarithmic coefficients"
            << nl << abort(FatalError);
    }

    intPolyType poly = this->integrate();

    return poly.evaluate(x2) - poly.evaluate(x1);
}


template<int PolySize>
typename Foam::Polynomial<PolySize>::intPolyType
Foam::Polynomial<PolySize>::integrate(const scalar intConstant)
{
    intPolyType newCoeffs;

    newCoeffs[0] = intConstant;
    forAll(*this, i)
    {
        newCoeffs[i + 1] = this->v_[i]/(i + 1);
    }

    return newCoeffs;
}


template<int PolySize>
typename Foam::Polynomial<PolySize>::polyType
Foam::Polynomial<PolySize>::integrateMinus1(const scalar intConstant)
{
    polyType newCoeffs;

    if (this->v_[0] > VSMALL)
    {
        newCoeffs.logActive() = true;
        newCoeffs.logCoeff() = this->v_[0];
    }

    newCoeffs[0] = intConstant;

    if (PolySize > 0)
    {
        for (label i=1; i<PolySize; i++)
        {
            newCoeffs[i] = this->v_[i]/i;
        }
    }

    return newCoeffs;
}


// ************************************************************************* //
