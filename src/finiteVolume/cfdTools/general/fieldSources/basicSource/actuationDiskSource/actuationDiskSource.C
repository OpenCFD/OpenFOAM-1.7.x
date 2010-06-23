/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "actuationDiskSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(actuationDiskSource, 0);
    addToRunTimeSelectionTable(basicSource, actuationDiskSource, dictionary);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::actuationDiskSource::checkData()
{
    if (magSqr(diskArea_) <= VSMALL)
    {
        FatalErrorIn("Foam::actuationDiskSource::checkData()")
           << "diskArea is approximately zero"
           << exit(FatalIOError);
    }
    if (Cp_ <= VSMALL || Ct_ <= VSMALL)
    {
        FatalErrorIn("Foam::actuationDiskSource::checkData()")
           << "Cp and Ct must be greater than zero"
           << exit(FatalIOError);
    }
    if (mag(diskDir_) < VSMALL)
    {
        FatalErrorIn("Foam::actuationDiskSource::checkData()")
           << "disk direction vector is approximately zero"
           << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::actuationDiskSource::actuationDiskSource
(
    const word& name,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    basicSource(name, dict, mesh),
    cellZoneID_(mesh.cellZones().findZoneID(this->cellSetName())),
    dict_(dict.subDict(typeName + "Coeffs")),
    diskDir_(dict_.lookup("diskDir")),
    Cp_(readScalar(dict_.lookup("Cp"))),
    Ct_(readScalar(dict_.lookup("Ct"))),
    diskArea_(readScalar(dict_.lookup("diskArea")))
{
    Info<< "    - creating actuation disk zone: "
        << this->name() << endl;

    bool foundZone = (cellZoneID_ != -1);

    reduce(foundZone, orOp<bool>());

    if (!foundZone && Pstream::master())
    {
        FatalErrorIn
        (
            "Foam::actuationDiskSource::actuationDiskSource"
            "(const word&, const dictionary&, const fvMesh&)"
        )   << "cannot find porous cellZone " << this->name()
            << exit(FatalError);
    }

    checkData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::actuationDiskSource::addSu(fvMatrix<vector>& UEqn)
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    bool compressible = false;
    if (UEqn.dimensions() == dimensionSet(1, 1, -2, 0, 0))
    {
        compressible = true;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = this->mesh().V();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();

    if (compressible)
    {
        addActuationDiskAxialInertialResistance
        (
            Usource,
            cells,
            V,
            this->mesh().lookupObject<volScalarField>("rho"),
            U
        );
    }
    else
    {
        addActuationDiskAxialInertialResistance
        (
            Usource,
            cells,
            V,
            geometricOneField(),
            U
        );
    }
}


void Foam::actuationDiskSource::writeData(Ostream& os) const
{

    os  << indent << token::BEGIN_BLOCK << incrIndent << nl;
    os.writeKeyword("name") << this->name() << token::END_STATEMENT << nl;

    if (dict_.found("note"))
    {
        os.writeKeyword("note") << string(dict_.lookup("note"))
            << token::END_STATEMENT << nl;
    }

    os << indent << "actuationDisk";

    dict_.write(os);

    os << decrIndent << indent << token::END_BLOCK << endl;
}


bool Foam::actuationDiskSource::read(const dictionary& dict)
{
    if (basicSource::read(dict))
    {
        const dictionary& sourceDict = dict.subDict(name());
        const dictionary& subDictCoeffs =
            sourceDict.subDict(typeName + "Coeffs");
        subDictCoeffs.readIfPresent("diskDir", diskDir_);
        subDictCoeffs.readIfPresent("Cp", Cp_);
        subDictCoeffs.readIfPresent("Ct", Ct_);
        subDictCoeffs.readIfPresent("diskArea", diskArea_);

        checkData();

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
