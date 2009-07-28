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

Application
    adiabaticFlameT

Description
    Calculates the adiabatic flame temperature for a given fuel over a
    range of unburnt temperatures and equivalence ratios.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OSspecific.H"

#include "specieThermo.H"
#include "janafThermo.H"
#include "perfectGas.H"

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
            << exit(FatalError);
    }

    dictionary control(controlFile);


    scalar T0(readScalar(control.lookup("T0")));
    word fuelName(control.lookup("fuel"));
    scalar n(readScalar(control.lookup("n")));
    scalar m(readScalar(control.lookup("m")));


    Info<< nl << "Reading Burcat data dictionary" << endl;

    fileName BurcatCpDataFileName(findEtcFile("thermoData/BurcatCpData"));

    // Construct control dictionary
    IFstream BurcatCpDataFile(BurcatCpDataFileName);

    // Check BurcatCpData stream is OK
    if (!BurcatCpDataFile.good())
    {
        FatalErrorIn(args.executable())
            << "Cannot read file " << BurcatCpDataFileName
            << exit(FatalError);
    }

    dictionary CpData(BurcatCpDataFile);


    scalar stoicO2 = n + m/4.0;
    scalar stoicN2 = (0.79/0.21)*(n + m/4.0);
    scalar stoicCO2 = n;
    scalar stoicH2O = m/2.0;

    thermo fuel
    (
        "fuel",
        thermo(CpData.lookup(fuelName))
    );

    thermo oxidant
    (
        "oxidant",
        stoicO2*thermo(CpData.lookup("O2"))
      + stoicN2*thermo(CpData.lookup("N2"))
    );

    dimensionedScalar stoichiometricAirFuelMassRatio
    (
        "stoichiometricAirFuelMassRatio",
        dimless,
        (oxidant.W()*oxidant.nMoles())/fuel.W()
    );

    Info<< "stoichiometricAirFuelMassRatio "
        << stoichiometricAirFuelMassRatio << ';' << endl;

    for (int i=0; i<300; i++)
    {
        scalar equiv = (i + 1)*0.01;
        scalar ft = 1/(1 + stoichiometricAirFuelMassRatio.value()/equiv);

        Info<< "phi = " << equiv << nl
            << "ft = " << ft << endl;

        scalar o2 = (1.0/equiv)*stoicO2;
        scalar n2 = (0.79/0.21)*o2;
        scalar fres = max(1.0 - 1.0/equiv, 0.0);
        scalar ores = max(1.0/equiv - 1.0, 0.0);
        scalar fburnt = 1.0 - fres;

        thermo fuel
        (
            "fuel",
            thermo(CpData.lookup(fuelName))
        );
        Info<< "fuel " << fuel << ';' << endl;

        thermo oxidant
        (
            "oxidant",
            o2*thermo(CpData.lookup("O2"))
          + n2*thermo(CpData.lookup("N2"))
        );
        Info<< "oxidant " << (1/oxidant.nMoles())*oxidant << ';' << endl;

        thermo reactants
        (
            "reactants",
            fuel + oxidant
        );
        Info<< "reactants " << (1/reactants.nMoles())*reactants << ';' << endl;

        thermo burntProducts
        (
            "burntProducts",
          + (n2 - (0.79/0.21)*ores*stoicO2)*thermo(CpData.lookup("N2"))
          + fburnt*stoicCO2*thermo(CpData.lookup("CO2"))
          + fburnt*stoicH2O*thermo(CpData.lookup("H2O"))
        );
        Info<< "burntProducts "
            << (1/burntProducts.nMoles())*burntProducts << ';' << endl;

        thermo products
        (
            "products",
            fres*fuel
          + n2*thermo(CpData.lookup("N2"))
          + fburnt*stoicCO2*thermo(CpData.lookup("CO2"))
          + fburnt*stoicH2O*thermo(CpData.lookup("H2O"))
          + ores*stoicO2*thermo(CpData.lookup("O2"))
        );

        Info<< "products " << (1/products.nMoles())*products << ';' << endl;

        scalar Tad = products.TH(reactants.H(T0), 1000.0);
        Info<< "Tad = " << Tad << nl << endl;
    }

    Info<< nl << "end" << endl;

    return 0;
}


// ************************************************************************* //
