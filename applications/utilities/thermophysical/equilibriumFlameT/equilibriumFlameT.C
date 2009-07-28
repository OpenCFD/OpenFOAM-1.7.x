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
    Calculates the equilibrium flame temperature for a given fuel and
    pressure for a range of unburnt gas temperatures and equivalence
    ratios; the effects of dissociation on O2, H2O and CO2 are included.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dictionary.H"
#include "IFstream.H"
#include "OSspecific.H"
#include "IOmanip.H"

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
            << abort(FatalError);
    }

    dictionary control(controlFile);


    scalar P(readScalar(control.lookup("P")));
    word fuel(control.lookup("fuel"));
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
            << abort(FatalError);
    }

    dictionary thermoData(BurcatCpDataFile);


    Info<< nl << "Reading Burcat data for relevant species" << nl << endl;

    // Reactants
    thermo FUEL(thermoData.lookup(fuel));
    thermo O2(thermoData.lookup("O2"));
    thermo N2(thermoData.lookup("N2"));

    // Products
    thermo CO2(thermoData.lookup("CO2"));
    thermo H2O(thermoData.lookup("H2O"));

    // Product fragments
    thermo CO(thermoData.lookup("CO"));
    thermo H2(thermoData.lookup("H2"));


    // Product dissociation reactions

    thermo CO2BreakUp
    (
        CO2 == CO + 0.5* O2
    );

    thermo H2OBreakUp
    (
        H2O == H2 + 0.5*O2
    );


    // Stoiciometric number of moles of species for one mole of fuel
    scalar stoicO2 = n + m/4.0;
    scalar stoicN2 = (0.79/0.21)*(n + m/4.0);
    scalar stoicCO2 = n;
    scalar stoicH2O = m/2.0;

    // Oxidant
    thermo oxidant
    (
        "oxidant",
        stoicO2*O2
      + stoicN2*N2
    );

    dimensionedScalar stoichiometricAirFuelMassRatio
    (
        "stoichiometricAirFuelMassRatio",
        dimless,
        (oxidant.W()*oxidant.nMoles())/FUEL.W()
    );

    Info<< "stoichiometricAirFuelMassRatio "
        << stoichiometricAirFuelMassRatio << ';' << endl;

    Info<< "Equilibrium flame temperature data ("
        << P/1e5 << " bar)" << nl << nl
        << setw(3) << "Phi"
        << setw(12) << "ft"
        << setw(7) << "T0"
        << setw(12) << "Tad"
        << setw(12) << "Teq"
        << setw(12) << "Terror"
        << setw(20) << "O2res (mole frac)" << nl
        << endl;


    // Loop over equivalence ratios
    for (int i=0; i<16; i++)
    {
        scalar equiv = 0.6 + i*0.05;
        scalar ft = 1/(1 + stoichiometricAirFuelMassRatio.value()/equiv);

    // Loop over initial temperatures
    for (int j=0; j<28; j++)
    {
        scalar T0 = 300.0 + j*100.0;

        // Number of moles of species for one mole of fuel
        scalar o2 = (1.0/equiv)*stoicO2;
        scalar n2 = (0.79/0.21)*o2;
        scalar fres = max(1.0 - 1.0/equiv, 0.0);
        scalar fburnt = 1.0 - fres;

        // Initial guess for number of moles of product species
        // ignoring product dissociation
        scalar oresInit = max(1.0/equiv - 1.0, 0.0)*stoicO2;
        scalar co2Init = fburnt*stoicCO2;
        scalar h2oInit = fburnt*stoicH2O;

        scalar ores = oresInit;
        scalar co2 = co2Init;
        scalar h2o = h2oInit;

        scalar co = 0.0;
        scalar h2 = 0.0;

        // Total number of moles in system
        scalar N = fres + n2 + co2 + h2o + ores;


        // Initial guess for adiabatic flame temperature
        scalar adiabaticFlameTemperature =
            T0
          + (fburnt/(1.0 + o2 + n2))/(1.0/(1.0 + (1.0 + 0.79/0.21)*stoicO2))
           *2000.0;

        scalar equilibriumFlameTemperature = adiabaticFlameTemperature;


        // Iteration loop for adiabatic flame temperature
        for (int j=0; j<20; j++)
        {

            if (j > 0)
            {
                co = co2*
                    min
                    (
                        CO2BreakUp.Kn(equilibriumFlameTemperature, P, N)
                       /::sqrt(max(ores, 0.001)),
                        1.0
                    );

                h2 = h2o*
                    min
                    (
                        H2OBreakUp.Kn(equilibriumFlameTemperature, P, N)
                       /::sqrt(max(ores, 0.001)),
                        1.0
                    );

                co2 = co2Init - co;
                h2o = h2oInit - h2;
                ores = oresInit + 0.5*co + 0.5*h2;
            }

            thermo reactants
            (
                FUEL + o2*O2 + n2*N2
            );

            thermo products
            (
                fres*FUEL + ores*O2 + n2*N2
              + co2*CO2 + h2o*H2O + co*CO + h2*H2
            );


            scalar equilibriumFlameTemperatureNew =
                products.TH(reactants.H(T0), adiabaticFlameTemperature);

            if (j==0)
            {
                adiabaticFlameTemperature = equilibriumFlameTemperatureNew;
            }
            else
            {
                equilibriumFlameTemperature = 0.5*
                (
                    equilibriumFlameTemperature
                  + equilibriumFlameTemperatureNew
                );
            }
        }

        Info<< setw(3) << equiv
            << setw(12) << ft
            << setw(7) << T0
            << setw(12) << adiabaticFlameTemperature
            << setw(12) << equilibriumFlameTemperature
            << setw(12)
            << adiabaticFlameTemperature - equilibriumFlameTemperature
            << setw(12) << ores/N
            << endl;
    }
    }

    Info<< nl << "end" << endl;

    return 0;
}


// ************************************************************************* //
