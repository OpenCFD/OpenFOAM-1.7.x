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

Application
    mixtureAdiabaticFlameT

Description
    Calculates the adiabatic flame temperature for a given mixture
    at a given temperature.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OSspecific.H"

#include "specieThermo.H"
#include "janafThermo.H"
#include "perfectGas.H"
#include "mixture.H"

using namespace Foam;

typedef specieThermo<janafThermo<perfectGas> > thermo;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.clear();
    argList::validArgs.append("controlFile");
    argList args(argc, argv);

    fileName controlFileName(args.additionalArgs()[0]);

    // Construct control dictionary
    IFstream controlFile(controlFileName);

    // Check controlFile stream is OK
    if (!controlFile.good())
    {
        FatalErrorIn(args.executable())
            << "Cannot read file " << controlFileName
            << abort(FatalError);
    }

    dictionary control(controlFile);


    scalar T0(readScalar(control.lookup("T0")));
    mixture rMix(control.lookup("reactants"));
    mixture pMix(control.lookup("products"));


    Info<< nl << "Reading Burcat data dictionary" << endl;

    fileName BurcatCpDataFileName(findEtcFile("thermoData/BurcatCpData"));

    // Construct control dictionary
    IFstream BurcatCpDataFile(BurcatCpDataFileName);

    // Check BurcatCpData stream is OK
    if (!BurcatCpDataFile.good())
    {
        FatalErrorIn(args.executable())
            << "Cannot read file " << BurcatCpDataFileName
            << abort(FatalError);
    }

    dictionary CpData(BurcatCpDataFile);


    thermo reactants
    (
        rMix[0].volFrac()*thermo(CpData.lookup(rMix[0].name()))
    );

    for (label i = 1; i < rMix.size(); i++)
    {
        reactants = reactants
            + rMix[i].volFrac()*thermo(CpData.lookup(rMix[i].name()));
    }


    thermo products
    (
        2*pMix[0].volFrac()*thermo(CpData.lookup(pMix[0].name()))
    );

    for (label i = 1; i < pMix.size(); i++)
    {
        products = products
            + 2*pMix[i].volFrac()*thermo(CpData.lookup(pMix[i].name()));
    }

    Info << "Adiabatic flame temperature of mixture " << rMix.name() << " = "
         << products.TH(reactants.H(T0), 1000.0) << " K" << endl;

    return 0;
}


// ************************************************************************* //
