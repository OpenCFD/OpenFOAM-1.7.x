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
    IFC (infinitely-fast chemistry) look-up table generator

Description
    Calculate the the infinitely-fast chemistry relationships in function of ft.
    for a given fuel.
    The output is given in moles.

    i.e. dictionary:

    fileName "SpeciesTable";


    fuel CH4(ANHARMONIC);
    n    1;
    m    4;


    fields
    (
        {
            name   ft;
            min    0.;
            max    1.;
            N      100;
        }
    );

    output
    (
        {
            name    CH4;
        }
        {
            name    CO2;
        }
        {
            name    H2O;
        }
);

\*---------------------------------------------------------------------------*/

#include "argList.H"

#include "IFstream.H"
#include "OFstream.H"

#include "specieThermo.H"
#include "janafThermo.H"
#include "perfectGas.H"

#include "IOdictionary.H"

#include "interpolationLookUpTable.H"

using namespace Foam;

typedef specieThermo<janafThermo<perfectGas> > thermo;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validArgs.clear();
    argList::validArgs.append("controlFile");
    argList args(argc, argv);

    fileName controlFileName(args.additionalArgs()[0]);

    enum varinput
    {
        fti,
        CH4i,
        CO2i,
        H2Oi
    };

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

    //  fuel-air mix
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

    word fuelName(control.lookup("fuel"));
    scalar n(readScalar(control.lookup("n")));
    scalar m(readScalar(control.lookup("m")));

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

    dimensionedScalar stoicRatio
    (
        "stoichiometricAirFuelMassRatio",
        dimless,
        (oxidant.W()*oxidant.nMoles())/(fuel.W()*fuel.nMoles())
    );


    // Open File for Look Up Table
    fileName LookUpTableFile(control.lookup("fileName"));

    OFstream controlFileOutput(LookUpTableFile);

    if (!controlFileOutput.good())
    {
        FatalErrorIn(args.executable())
            << "Cannot open file " << LookUpTableFile
            << exit(FatalError);
    }

    // Create Look Up Table
    interpolationLookUpTable<scalar> LookUpTable(control);

    const List<label>& dim = LookUpTable.dim();

    const List<scalar>& min = LookUpTable.min();

    const List<scalar>& delta = LookUpTable.delta();

    label count = 0;

    for (label i=0; i <= dim[fti]; i++)
    {
        scalar ft = Foam::min(scalar(i)*delta[fti] + min[fti] + 0.001, 0.999);

        scalar equiv = Foam::pow(((1.0 / ft) - 1.0), -1.0)*stoicRatio.value();

        scalar o2 = (1.0/equiv)*stoicO2;
        scalar n2 = (0.79/0.21)*o2;
        scalar fres = max(1.0 - 1.0/equiv, 0.0);
        scalar ores = max(1.0/equiv - 1.0, 0.0);
        scalar fburnt = 1.0 - fres;


        thermo fuel
        (
            "fuel",
            fres*thermo(CpData.lookup(fuelName))
        );

        thermo N2
        (
            "N2",
            n2*thermo(CpData.lookup("N2"))
        );

        thermo O2
        (
            "O2",
            ores*thermo(CpData.lookup("O2"))
        );

        thermo CO2
        (
            "CO2",
            fburnt*stoicCO2*thermo(CpData.lookup("CO2"))
        );

        thermo H2O
        (
            "H2O",
            fburnt*stoicH2O*thermo(CpData.lookup("H2O"))
        );


        scalar ToTalMoles = fuel.nMoles() + CO2.nMoles() + H2O.nMoles() +
                N2.nMoles() + O2.nMoles();

        LookUpTable[fti][count] = ft;
        LookUpTable[CH4i][count] = fuel.nMoles()/ToTalMoles;
        LookUpTable[CO2i][count] = CO2.nMoles()/ToTalMoles;
        LookUpTable[H2Oi][count] = H2O.nMoles()/ToTalMoles;
        count++;
    }

    IOobject::writeBanner(controlFileOutput);
    controlFileOutput << "\n" << nl;
    controlFileOutput.writeKeyword("fields");
    controlFileOutput << LookUpTable.entries() << token::END_STATEMENT << nl;

    controlFileOutput.writeKeyword("output");
    controlFileOutput << LookUpTable.output() << token::END_STATEMENT << nl;

    if (LookUpTable.size() == 0)
    {
        FatalErrorIn
        (
            "Foam::IFCLookUpTableGen"
        )   << "table is empty" << nl
            << exit(FatalError);
    }

    controlFileOutput.writeKeyword("values");
    controlFileOutput << LookUpTable << token::END_STATEMENT << nl;

    return(0);
}


// ************************************************************************* //
