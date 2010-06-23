/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "basicSourceList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicSourceList::basicSourceList
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    PtrList<basicSource>(),
    mesh_(mesh)
{
    label count = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        // safety:
        if (iter().isDict())
        {
            count ++;
        }
    }

    this->setSize(count);
    label i = 0;
    forAllConstIter(dictionary, dict, iter)
    {
        const word& name = iter().keyword();
        const dictionary& dict = iter().dict();

        this->set
        (
            i++,
            basicSource::New(name, dict, mesh)
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::basicSourceList::addSu(fvMatrix<scalar>& Eqn)
{
    forAll(*this, i)
    {
        if (this->operator[](i).isActive())
        {
            this->operator[](i).addSu(Eqn);
        }
    }
}


void Foam::basicSourceList::addSu(fvMatrix<vector>& Eqn)
{

    forAll(*this, i)
    {
        if (this->operator[](i).isActive())
        {
            this->operator[](i).addSu(Eqn);
        }
    }
}


void Foam::basicSourceList::addExplicitSources()
{

    forAll(*this, i)
    {
        if (this->operator[](i).isActive())
        {
            this->operator[](i).addExplicitSources();
        }
    }
}


void Foam::basicSourceList::addSu
(
    DimensionedField<scalar, volMesh>& field
)
{
    forAll(*this, i)
    {
        if (this->operator[](i).isActive())
        {
            this->operator[](i).addSu(field);
        }
    }
}


void Foam::basicSourceList::addSu
(
    DimensionedField<vector, volMesh>& field
)
{
    forAll(*this, i)
    {
        if (this->operator[](i).isActive())
        {
            this->operator[](i).addSu(field);
        }
    }
}


bool Foam::basicSourceList::read(const dictionary& dict)
{
    forAll(*this, i)
    {
        this->operator[](i).read(dict);
    }
    return true;
}


bool Foam::basicSourceList::writeData(Ostream& os) const
{
    // Write size of list
    os << nl << this->size();

    // Write beginning of contents
    os << nl << token::BEGIN_LIST;

    // Write list contents
    forAll(*this, i)
    {
        os << nl;
        this->operator[](i).writeData(os);
    }

    // Write end of contents
    os << token::END_LIST << token::END_STATEMENT << nl;

    // Check state of IOstream
     return os.good();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const basicSourceList& sources
)
{
    sources.writeData(os);
    return os;
}


// ************************************************************************* //
