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

#include "error.H"

#include "blobsSheetAtomization.H"
#include "addToRunTimeSelectionTable.H"
#include "basicMultiComponentMixture.H"

#include "RosinRammler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(blobsSheetAtomization, 0);

addToRunTimeSelectionTable
(
    atomizationModel,
    blobsSheetAtomization,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
blobsSheetAtomization::blobsSheetAtomization
(
    const dictionary& dict,
    spray& sm
)
:
    atomizationModel(dict, sm),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    B_(readScalar(coeffsDict_.lookup("B"))),
    angle_(readScalar(coeffsDict_.lookup("angle"))),
    rndGen_(sm.rndGen())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

blobsSheetAtomization::~blobsSheetAtomization()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void blobsSheetAtomization::atomizeParcel
(
    parcel& p,
    const scalar deltaT,
    const vector& vel,
    const liquidMixture& fuels
) const
{

    const PtrList<volScalarField>& Y = spray_.composition().Y();

    label Ns = Y.size();
    label cellI = p.cell();
    scalar pressure = spray_.p()[cellI];
    scalar temperature = spray_.T()[cellI];
    scalar Taverage = p.T() + (temperature - p.T())/3.0;

    scalar Winv = 0.0;
    for(label i=0; i<Ns; i++)
    {
        Winv += Y[i][cellI]/spray_.gasProperties()[i].W();
    }
    scalar R = specie::RR*Winv;

    // ideal gas law to evaluate density
    scalar rhoAverage = pressure/R/Taverage;
    scalar sigma = fuels.sigma(pressure, p.T(), p.X());

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //     The We and Re numbers are to be evaluated using the 1/3 rule.
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    scalar rhoFuel = fuels.rho(1.0e+5, p.T(), p.X());

    scalar U = mag(p.Urel(vel));

    const injectorType& it =
        spray_.injectors()[label(p.injector())].properties();

    vector itPosition(vector::zero);
    label nHoles = it.nHoles();
    if (nHoles > 1)
    {
        for(label i=0; i<nHoles;i++)
        {
            itPosition += it.position(i);
        }
        itPosition /= nHoles;
    }
    else
    {
        itPosition = it.position(0);
    }
//    const vector itPosition = it.position();


    scalar lBU = B_ * sqrt
    (
        rhoFuel * sigma * p.d() * cos(angle_*mathematicalConstant::pi/360.0)
      / sqr(rhoAverage*U)
    );

    scalar pWalk = mag(p.position() - itPosition);

    if(pWalk > lBU && p.liquidCore() == 1.0)
    {
        p.liquidCore() = 0.0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
